"""Structure geometry for proximity-based design: ligand burial, exit vectors, ternary poses.

What this does, concretely:
  1. fetches a PDB entry (mmCIF, cached),
  2. measures how buried the bound ligand is (Shrake-Rupley SASA, per atom),
  3. ranks the ligand's solvent-exposed atoms as linker attachment points ("exit vectors"),
  4. for a pair of ligand-bound structures, places the second protein so the two exit
     vectors face each other at a given linker length, spins it about that axis and
     reports clashes and interface contacts for every pose.

What this is NOT: ternary-complex prediction. Step 4 is rigid-body geometric triage of a
single crystal conformation. It answers "could these two ligands point at each other with
this much linker without the proteins overlapping", not "will this PROTAC work". Real
ternary complexes involve induced fit, linker chemistry and cooperativity that no rigid
superposition captures. Treat a clash-free pose as a hypothesis to model properly, and a
badly clashing one as a reason to look at a different exit vector or a longer linker.
"""
import gzip
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from Bio.PDB import MMCIFParser, PDBIO, Select
from Bio.PDB.kdtrees import KDTree
from Bio.PDB.SASA import ShrakeRupley

from .mechanisms import CACHE
from .structures import ADDITIVES, ligand_kind

RCSB_FILE = "https://files.rcsb.org/download/{pdb_id}.cif.gz"
CLASH = 2.5        # heavy-atom distance below this = steric clash
CONTACT = 4.5      # heavy-atom distance below this = interface contact
LINKERS = (5.0, 8.0, 11.0, 14.0, 17.0, 20.0)   # attachment-point separations to scan, A


# ---------- loading ----------

def structure_path(pdb_id, refresh=False):
    d = Path(CACHE) / "structures"
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"{pdb_id.upper()}.cif"
    if refresh or not path.exists() or path.stat().st_size == 0:
        r = requests.get(RCSB_FILE.format(pdb_id=pdb_id.upper()), timeout=180)
        r.raise_for_status()
        path.write_bytes(gzip.decompress(r.content))
    return path


def load(pdb_id, refresh=False):
    """First model of a PDB entry, hydrogens and waters left as they are."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = MMCIFParser(QUIET=True).get_structure(pdb_id.upper(), str(structure_path(pdb_id, refresh)))
    return next(s.get_models())


def _heavy(atoms):
    return [a for a in atoms if a.element not in ("H", "D")]


def ligands(model, min_atoms=6):
    """Non-additive hetero residues in the model, biggest first."""
    out = []
    for chain in model:
        for res in chain:
            het = res.id[0]
            if not het.startswith("H_"):
                continue
            name = res.get_resname().strip()
            if name in ADDITIVES:
                continue
            heavy = _heavy(res.get_atoms())
            if len(heavy) < min_atoms:
                continue
            mw = sum(_MASS.get(a.element, 12.0) for a in heavy)
            out.append({"resname": name, "chain": chain.id, "resseq": res.id[1],
                        "n_atoms": len(heavy), "approx_mw": round(mw, 1),
                        "kind": ligand_kind(name, mw), "residue": res})
    return sorted(out, key=lambda d: -d["n_atoms"])


def pick_ligand(model, ligand_id=None, chain=None):
    """The ligand to work with: named one if given, else the largest non-additive one."""
    cand = ligands(model)
    if ligand_id:
        cand = [c for c in cand if c["resname"] == ligand_id.upper()]
    if chain:
        cand = [c for c in cand if c["chain"] == chain]
    if not cand:
        raise ValueError(f"no suitable ligand found (asked for {ligand_id or 'any'})")
    return cand[0]


_MASS = {"C": 12.011, "N": 14.007, "O": 15.999, "S": 32.06, "P": 30.974, "F": 18.998,
         "CL": 35.45, "BR": 79.904, "I": 126.904, "B": 10.81, "SI": 28.085}


# ---------- burial ----------

def ligand_sasa(model, lig, probe_radius=1.4, n_points=200):
    """Per-atom and total ligand SASA, in the complex and for the free ligand.

    buried_fraction = 1 - SASA(in complex) / SASA(ligand alone). 1.0 = fully enclosed.
    """
    res = lig["residue"]
    sr = ShrakeRupley(probe_radius=probe_radius, n_points=n_points)

    sr.compute(model, level="A")
    in_complex = {a.get_id(): float(a.sasa) for a in res.get_atoms()}

    free_struct = res.copy()                       # detached copy: ligand on its own
    free_struct.detach_parent()
    sr.compute(free_struct, level="A")
    free = {a.get_id(): float(a.sasa) for a in free_struct.get_atoms()}

    rows = []
    for a in _heavy(res.get_atoms()):
        aid = a.get_id()
        f = free.get(aid, 0.0)
        c = in_complex.get(aid, 0.0)
        rows.append({"atom": aid, "element": a.element,
                     "sasa_complex": round(c, 2), "sasa_free": round(f, 2),
                     "rel_exposure": round(c / f, 3) if f > 1e-6 else 0.0,
                     "coord": a.get_coord()})
    tot_c = sum(r["sasa_complex"] for r in rows)
    tot_f = sum(r["sasa_free"] for r in rows)
    return {"atoms": rows, "sasa_complex": round(tot_c, 1), "sasa_free": round(tot_f, 1),
            "buried_fraction": round(1 - tot_c / tot_f, 3) if tot_f > 1e-6 else None,
            "n_exposed_atoms": sum(1 for r in rows if r["rel_exposure"] >= 0.2)}


# ---------- exit vectors ----------

def _unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-9 else np.array([0.0, 0.0, 1.0])


def _protein_atoms(model, exclude_res, keep_chains=None):
    out = []
    for chain in model:
        if keep_chains is not None and chain.id not in keep_chains:
            continue
        for res in chain:
            if res is exclude_res or res.id[0] not in (" ",):
                continue
            out += _heavy(res.get_atoms())
    return out


def exit_vectors(model, lig, sasa=None, min_rel=0.15, cone_deg=35.0, probe_len=12.0, top=8,
                 keep_chains=None, fallback=True):
    """Rank the ligand's exposed atoms as linker attachment points.

    direction: away from the local protein wall (atom - centre of protein atoms within 8 A).
    openness:  fraction of a cone (half-angle cone_deg, length probe_len) along that
               direction that is free of protein atoms -> 1.0 = nothing in the way.

    A ligand in a fully enclosed pocket has no atom above `min_rel` (BRAF's WQN in 8VSO is
    99.6% buried). With fallback=True the most exposed atoms are returned anyway, flagged
    `buried=True`: that is a real answer for proximity design - this structure offers no
    exit vector, and the pocket would have to open for a linker to leave.
    """
    sasa = sasa or ligand_sasa(model, lig)
    prot = _protein_atoms(model, lig["residue"], keep_chains)
    if not prot:
        raise ValueError("structure has no protein atoms around the ligand")
    pcoords = np.array([a.get_coord() for a in prot], dtype="d")
    kdt = KDTree(pcoords, 10)
    lig_centre = np.mean([r["coord"] for r in sasa["atoms"]], axis=0)
    cos_cut = np.cos(np.radians(cone_deg))

    chosen = [r for r in sasa["atoms"] if r["rel_exposure"] >= min_rel]
    buried = not chosen
    if buried and fallback:
        chosen = sorted(sasa["atoms"], key=lambda r: -r["rel_exposure"])[:top]
    rows = []
    for r in chosen:
        pos = np.asarray(r["coord"], dtype="d")
        near = [p.index for p in kdt.search(pos, 8.0)]
        wall = pcoords[near].mean(axis=0) if near else lig_centre
        direction = _unit(pos - wall)
        if np.dot(direction, _unit(pos - lig_centre)) < 0:   # keep it pointing outward
            direction = _unit(pos - lig_centre)

        blocked = 0
        for p in kdt.search(pos, probe_len):
            d = pcoords[p.index] - pos
            dist = np.linalg.norm(d)
            if dist < 1e-6:
                continue
            if np.dot(d / dist, direction) >= cos_cut:
                blocked += 1
        rows.append({"atom": r["atom"], "element": r["element"],
                     "rel_exposure": r["rel_exposure"], "sasa_complex": r["sasa_complex"],
                     "blocked_atoms": blocked, "openness": round(1.0 / (1.0 + blocked), 3),
                     "buried": buried, "coord": pos, "direction": direction})
    rows.sort(key=lambda d: (-d["openness"], -d["rel_exposure"]))
    return rows[:top]


def exit_vector_table(rows, pdb_id="", ligand=""):
    return pd.DataFrame([{"pdb_id": pdb_id, "ligand": ligand, "atom": r["atom"],
                          "element": r["element"], "rel_exposure": r["rel_exposure"],
                          "sasa_complex": r["sasa_complex"], "blocked_atoms": r["blocked_atoms"],
                          "openness": r["openness"],
                          "no_open_vector": bool(r.get("buried")),
                          "x": round(float(r["coord"][0]), 2),
                          "y": round(float(r["coord"][1]), 2),
                          "z": round(float(r["coord"][2]), 2)} for r in rows])


# ---------- rigid-body ternary placement ----------

def _rot_between(a, b):
    """Rotation matrix taking unit vector a onto unit vector b (Rodrigues)."""
    a, b = _unit(np.asarray(a, "d")), _unit(np.asarray(b, "d"))
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-9:
        if c > 0:
            return np.eye(3)
        axis = _unit(np.cross(a, [1.0, 0.0, 0.0]) if abs(a[0]) < 0.9
                     else np.cross(a, [0.0, 1.0, 0.0]))
        return _rot_axis(axis, np.pi)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]], dtype="d")
    return np.eye(3) + vx + vx @ vx * (1.0 / (1.0 + c))


def _rot_axis(axis, angle):
    axis = _unit(np.asarray(axis, "d"))
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]], "d")
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)


def _model_coords(model, keep_chains=None, include_het=False):
    atoms = []
    for chain in model:
        if keep_chains is not None and chain.id not in keep_chains:
            continue
        for res in chain:
            if res.id[0] != " " and not include_het:
                continue
            atoms += _heavy(res.get_atoms())
    return atoms, np.array([a.get_coord() for a in atoms], dtype="d")


def chains_near_ligand(model, lig, cutoff=6.0, min_frac=0.15):
    """Protein chains that really bind the ligand, as {chain_id}.

    PDB files often hold several copies of a complex, and a neighbouring copy may brush
    the ligand through crystal packing. Chains are ranked by how many atom pairs fall
    within `cutoff`; a chain is kept only if it reaches `min_frac` of the top chain
    (e.g. in 4W9H the real chain has 690 such pairs and a packing neighbour 13).
    Pass an explicit chain list to override.
    """
    lig_xyz = np.array([a.get_coord() for a in _heavy(lig["residue"].get_atoms())], dtype="d")
    counts = {}
    for chain in model:
        xyz = [a.get_coord() for res in chain if res.id[0] == " " for a in _heavy(res.get_atoms())]
        if not xyz:
            continue
        d = np.linalg.norm(np.asarray(xyz, "d")[:, None, :] - lig_xyz[None, :, :], axis=-1)
        n = int((d < cutoff).sum())
        if n:
            counts[chain.id] = n
    if not counts:
        return set()
    top = max(counts.values())
    return {c for c, n in counts.items() if n >= min_frac * top}


def _count_pairs(ca, cb, centre_a, radius_a, clash, contact):
    """Heavy-atom clash and contact counts between a fixed set A and a posed set B.

    Only B atoms that could possibly reach A (inside A's bounding sphere + contact) and
    then only A atoms near those are compared, so each pose costs a small matrix.
    """
    near_b = cb[np.linalg.norm(cb - centre_a, axis=1) <= radius_a + contact]
    if not len(near_b):
        return 0, 0
    centre_b = near_b.mean(axis=0)
    radius_b = float(np.linalg.norm(near_b - centre_b, axis=1).max())
    near_a = ca[np.linalg.norm(ca - centre_b, axis=1) <= radius_b + contact]
    if not len(near_a):
        return 0, 0
    n_clash = n_contact = 0
    for i in range(0, len(near_b), 256):
        d = np.linalg.norm(near_b[i:i + 256, None, :] - near_a[None, :, :], axis=-1)
        n_contact += int((d < contact).sum())
        n_clash += int((d < clash).sum())
    return n_clash, n_contact


def _apply(xyz, rot):
    """Rotate coordinates by `rot`.

    numpy 2.0 on Apple's Accelerate BLAS raises spurious divide/overflow/invalid warnings
    for perfectly ordinary float64 matmuls (even `a @ np.eye(3)`), so they are silenced
    here; the arithmetic itself is unaffected.
    """
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        return np.asarray(xyz) @ rot.T


def _perp(v):
    """Some unit vector perpendicular to v."""
    ref = np.array([1.0, 0.0, 0.0]) if abs(v[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    return _unit(np.cross(v, ref))


def _axis_directions(dir_a, tilts, n_azimuth):
    """Directions for the A-anchor -> B-anchor displacement: dir_a, then tilted rings.

    A linker is flexible and curved, so the partner does not have to sit exactly along
    the exit vector. Each tilt ring relaxes that by a few degrees.
    """
    out = [(_unit(dir_a), 0.0, 0.0)]
    p = _perp(dir_a)
    for tilt in tilts:
        if tilt <= 0:
            continue
        tilted = _rot_axis(p, np.radians(tilt)) @ dir_a
        for j in range(n_azimuth):
            phi = 2 * np.pi * j / n_azimuth
            out.append((_unit(_rot_axis(dir_a, phi) @ tilted), float(tilt),
                        round(np.degrees(phi), 1)))
    return out


def ternary_scan(model_a, lig_a, vec_a, model_b, lig_b, vec_b,
                 linkers=LINKERS, n_rot=36, clash=CLASH, contact=CONTACT,
                 chains_a=None, chains_b=None, tilts=(0,), n_azimuth=6):
    """Place B so its ligand exit vector faces A's, and scan the free parameters.

    Sampled: linker length (anchor separation), the tilt of the A->B displacement away
    from A's exit vector, the azimuth of that tilt, and the spin of B about the
    displacement axis. Returns one row per pose with clash and contact counts.
    """
    atoms_a, ca = _model_coords(model_a, chains_a)
    atoms_b, cb = _model_coords(model_b, chains_b)
    if not len(ca) or not len(cb):
        raise ValueError("one of the structures has no protein atoms")
    centre_a = ca.mean(axis=0)
    radius_a = float(np.linalg.norm(ca - centre_a, axis=1).max())

    anchor_a = np.asarray(vec_a["coord"], "d")
    dir_a = _unit(np.asarray(vec_a["direction"], "d"))
    anchor_b = np.asarray(vec_b["coord"], "d")
    dir_b = _unit(np.asarray(vec_b["direction"], "d"))

    rows = []
    for axis, tilt, azim in _axis_directions(dir_a, tilts, n_azimuth):
        Ra = _rot_between(dir_b, -axis)            # B's vector points back along the axis
        base_axis = _apply(cb - anchor_b, Ra)
        for L in linkers:
            target = anchor_a + axis * L
            for i in range(n_rot):
                ang = 2 * np.pi * i / n_rot
                spin = _rot_axis(axis, ang)
                coords = _apply(base_axis, spin) + target
                n_clash, n_contact = _count_pairs(ca, coords, centre_a, radius_a, clash, contact)
                rows.append({"linker_A": L, "tilt_deg": tilt, "azimuth_deg": azim,
                             "spin_deg": round(np.degrees(ang), 1),
                             "clashes": n_clash, "contacts": n_contact,
                             "_R": spin @ Ra, "_t": target, "_anchor_b": anchor_b,
                             "_centroid": coords.mean(axis=0)})
    rows.sort(key=lambda d: (d["clashes"], -d["contacts"]))
    return rows


def cluster_poses(rows, cutoff=8.0, max_clusters=10, clash_tol=0):
    """Distinct clash-free arrangements, best first.

    Thousands of poses differing by a few degrees are not thousands of hypotheses. Poses
    are taken best-first (fewest clashes, then most contacts) and a pose is kept only if
    the partner's centre sits more than `cutoff` from every pose already kept.
    """
    ok = [r for r in rows if r["clashes"] <= clash_tol]
    ok.sort(key=lambda d: (d["clashes"], -d["contacts"]))
    reps = []
    for r in ok:
        c = r["_centroid"]
        if all(np.linalg.norm(c - k["_centroid"]) > cutoff for k in reps):
            reps.append(r)
            if len(reps) >= max_clusters:
                break
    return reps


def scan_table(rows):
    return pd.DataFrame([{k: v for k, v in r.items() if not k.startswith("_")} for r in rows])


def write_pose(model_a, model_b, pose, path, chains_a=None, chains_b=None,
               chain_offset="abcdefghijklmnopqrstuvwxyz"):
    """Write A plus the transformed B as one PDB file, for PyMOL/ChimeraX."""
    import copy

    b = copy.deepcopy(model_b)
    for chain in list(b):
        if chains_b is not None and chain.id not in chains_b:
            b.detach_child(chain.id)
    R, t, anchor = pose["_R"], pose["_t"], pose["_anchor_b"]
    for atom in b.get_atoms():
        atom.set_coord(_apply(atom.get_coord() - anchor, R) + t)

    merged = copy.deepcopy(model_a)
    for chain in list(merged):
        if chains_a is not None and chain.id not in chains_a:
            merged.detach_child(chain.id)
    used = {c.id for c in merged}
    free = [c for c in chain_offset if c.upper() not in {u.upper() for u in used}]
    for i, chain in enumerate(list(b)):
        chain.detach_parent()
        chain.id = free[i] if i < len(free) else f"{i}"
        merged.add(chain)

    io = PDBIO()
    io.set_structure(merged)
    io.save(str(path), select=_NoWater())
    return path


class _NoWater(Select):
    def accept_residue(self, residue):
        return residue.get_resname().strip() not in ("HOH", "DOD")


# ---------- one call per pair ----------

def analyse_pair(pdb_a, pdb_b, ligand_a=None, ligand_b=None, n_vectors=3,
                 linkers=LINKERS, n_rot=36, out_dir=None, tag="",
                 chains_a=None, chains_b=None, tilts=(0, 20, 40), n_azimuth=6,
                 cluster_cutoff=8.0, max_poses=5):
    """Burial + exit vectors for both structures, then a pose scan over vector pairs.

    Returns (summary dict, vectors DataFrame, poses DataFrame).
    """
    ma, mb = load(pdb_a), load(pdb_b)
    la, lb = pick_ligand(ma, ligand_a), pick_ligand(mb, ligand_b)
    ka = set(chains_a) if chains_a else chains_near_ligand(ma, la)
    kb = set(chains_b) if chains_b else chains_near_ligand(mb, lb)
    sa, sb = ligand_sasa(ma, la), ligand_sasa(mb, lb)
    va = exit_vectors(ma, la, sa, top=n_vectors, keep_chains=ka)
    vb = exit_vectors(mb, lb, sb, top=n_vectors, keep_chains=kb)
    if not va or not vb:
        raise ValueError("ligand has no atoms to attach a linker to")
    warn = []
    if va[0].get("buried"):
        warn.append(f"{pdb_a.upper()}/{la['resname']} is {sa['buried_fraction']:.0%} buried")
    if vb[0].get("buried"):
        warn.append(f"{pdb_b.upper()}/{lb['resname']} is {sb['buried_fraction']:.0%} buried")

    poses = []
    for ia, ra in enumerate(va):
        for ib, rb in enumerate(vb):
            for row in ternary_scan(ma, la, ra, mb, lb, rb, linkers, n_rot,
                                    chains_a=ka, chains_b=kb, tilts=tilts,
                                    n_azimuth=n_azimuth):
                row.update({"atom_a": ra["atom"], "atom_b": rb["atom"],
                            "openness_a": ra["openness"], "openness_b": rb["openness"]})
                poses.append(row)
    poses.sort(key=lambda d: (d["clashes"], -d["contacts"]))
    reps = cluster_poses(poses, cutoff=cluster_cutoff, max_clusters=max_poses)

    n_free = int(sum(1 for p in poses if p["clashes"] == 0))
    best = reps[0] if reps else poses[0]
    files = []
    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        stem = tag or f"{pdb_a}_{pdb_b}"
        for i, r in enumerate(reps or [best], start=1):
            name = f"pose_{stem}_{i}.pdb"
            write_pose(ma, mb, r, out_dir / name, chains_a=ka, chains_b=kb)
            files.append(name)

    cols = ["atom_a", "atom_b", "linker_A", "tilt_deg", "azimuth_deg", "spin_deg",
            "clashes", "contacts", "openness_a", "openness_b"]
    pdf = pd.DataFrame([{k: p[k] for k in cols} for p in poses])
    summary = {
        "pdb_a": pdb_a.upper(), "ligand_a": la["resname"],
        "chains_a": ",".join(sorted(ka)), "pdb_b": pdb_b.upper(), "ligand_b": lb["resname"],
        "chains_b": ",".join(sorted(kb)),
        "buried_fraction_a": sa["buried_fraction"], "buried_fraction_b": sb["buried_fraction"],
        "sasa_complex_a": sa["sasa_complex"], "sasa_complex_b": sb["sasa_complex"],
        "n_exposed_atoms_a": sa["n_exposed_atoms"], "n_exposed_atoms_b": sb["n_exposed_atoms"],
        "top_atom_a": best["atom_a"], "top_atom_b": best["atom_b"],
        "top_linker_A": best["linker_A"], "top_clashes": best["clashes"],
        "top_contacts": best["contacts"],
        "clash_free_poses": n_free, "n_poses": len(poses),
        "clash_free_pct": round(100.0 * n_free / len(poses), 1) if poses else 0.0,
        "distinct_arrangements": len(reps),
        "linkers_with_clash_free": ", ".join(
            f"{L:g}" for L in sorted({p["linker_A"] for p in poses if p["clashes"] == 0})),
        "pose_files": ", ".join(files),
        "warning": ("; ".join(warn) + " - no solvent-exposed exit vector in this structure")
                   if warn else "",
    }
    vecs = pd.concat([exit_vector_table(va, pdb_a.upper(), la["resname"]),
                      exit_vector_table(vb, pdb_b.upper(), lb["resname"])], ignore_index=True)
    return summary, vecs, pdf
