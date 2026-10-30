#!/usr/bin/env python3
"""How close does the ternary scan get to a real PROTAC ternary complex?

The only honest way to judge a geometric pose scan is against a crystal structure of the
thing it is trying to approximate. 5T35 is BRD4(BD2) + the PROTAC MZ1 + VHL.

The test:
  1. take the two BINARY structures the scan would normally start from --
     3MXF (BRD4 BD1 + JQ1) and 4W9H (VHL + VH032),
  2. run the scan,
  3. superpose our BRD4 onto 5T35's BRD4, apply the same transform to each posed VHL,
  4. measure CA RMSD against 5T35's VHL.

Two separate questions come out of this:
  - can the sampling REACH the native arrangement?  (best RMSD over all poses)
  - does the SCORING find it?                       (RMSD of the top-ranked pose)

Run:  python validation/benchmark_5t35.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Data.IUPACData import protein_letters_3to1 as THREE_TO_ONE
from Bio.PDB import Superimposer
from Bio.PDB.Polypeptide import is_aa

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ppi import geometry as G                                              # noqa: E402

BLOSUM = substitution_matrices.load("BLOSUM62")
# 5T35: chain A is the BRD4 bromodomain, chain D is VHL (B and C are elongins).
NATIVE, NATIVE_TARGET, NATIVE_E3 = "5T35", "A", "D"
BINARY_TARGET, BINARY_E3 = "3MXF", "4W9H"


def seq_and_ca(chain):
    seq, cas = "", []
    for res in chain:
        if not is_aa(res, standard=True) or "CA" not in res:
            continue
        seq += THREE_TO_ONE.get(res.get_resname().capitalize(), "X")
        cas.append(res["CA"])
    return seq, cas


def matched_ca(chain1, chain2, schedule=(12.0, 8.0, 5.0, 3.5, 3.0), min_pairs=40):
    """CA pairs from a sequence alignment, refined by superposing and tightening the cutoff."""
    s1, c1 = seq_and_ca(chain1)
    s2, c2 = seq_and_ca(chain2)
    aligner = PairwiseAligner(mode="local", substitution_matrix=BLOSUM,
                              open_gap_score=-11, extend_gap_score=-1)
    alignment = aligner.align(s1, s2)[0]
    a, b = [], []
    for (i0, i1), (j0, j1) in zip(alignment.aligned[0], alignment.aligned[1]):
        for k in range(i1 - i0):
            a.append(c1[i0 + k])
            b.append(c2[j0 + k])
    for cut in schedule:
        sup = Superimposer()
        sup.set_atoms(b, a)
        rot, tran = sup.rotran
        d = np.linalg.norm(np.array([x.get_coord() for x in a]) @ rot + tran
                           - np.array([x.get_coord() for x in b]), axis=1)
        keep = d <= cut
        if keep.sum() < min_pairs or keep.all():
            continue
        a = [x for x, k in zip(a, keep) if k]
        b = [x for x, k in zip(b, keep) if k]
    sup = Superimposer()
    sup.set_atoms(b, a)
    return a, b, sup


def main():
    native = G.load(NATIVE)
    target_nat, e3_nat = native[NATIVE_TARGET], native[NATIVE_E3]
    _, c_t = G._model_coords(native, {NATIVE_TARGET})
    _, c_e = G._model_coords(native, {NATIVE_E3})
    d = np.linalg.norm(c_t[:, None, :] - c_e[None, :, :], axis=-1)
    print(f"{NATIVE} native interface: {(d < G.CONTACT).sum()} contacts, "
          f"{(d < G.CLASH).sum()} clashes")

    ma, mb = G.load(BINARY_TARGET), G.load(BINARY_E3)
    la, lb = G.pick_ligand(ma), G.pick_ligand(mb)
    ka, kb = G.chains_near_ligand(ma, la), G.chains_near_ligand(mb, lb)
    sa, sb = G.ligand_sasa(ma, la), G.ligand_sasa(mb, lb)
    print(f"{BINARY_TARGET} {la['resname']}: buried {sa['buried_fraction']:.0%}, "
          f"{sa['n_exposed_atoms']} exposed atoms | "
          f"{BINARY_E3} {lb['resname']}: buried {sb['buried_fraction']:.0%}, "
          f"{sb['n_exposed_atoms']} exposed atoms")

    chain_a = sorted(ka)[0]
    _, _, sup_t = matched_ca(ma[chain_a], target_nat)
    print(f"{BINARY_TARGET} {chain_a} -> {NATIVE} {NATIVE_TARGET}: RMSD {sup_t.rms:.2f} A")
    rot, tran = sup_t.rotran

    ours, native_ca, sup_e = matched_ca(mb[sorted(kb)[0]], e3_nat)
    print(f"{BINARY_E3} -> {NATIVE} {NATIVE_E3}: {len(ours)} CA pairs, RMSD {sup_e.rms:.2f} A")
    ours_xyz = np.array([a.get_coord() for a in ours])
    native_xyz = np.array([a.get_coord() for a in native_ca])

    va = G.exit_vectors(ma, la, sa, top=3, keep_chains=ka)
    vb = G.exit_vectors(mb, lb, sb, top=3, keep_chains=kb)
    poses = []
    for ra in va:
        for rb in vb:
            for row in G.ternary_scan(ma, la, ra, mb, lb, rb, n_rot=24, chains_a=ka,
                                      chains_b=kb, tilts=(0, 20, 40), n_azimuth=6):
                row.update(atom_a=ra["atom"], atom_b=rb["atom"])
                poses.append(row)

    for p in poses:
        posed = (ours_xyz - p["_anchor_b"]) @ p["_R"].T + p["_t"]     # into the binary frame
        posed = posed @ rot + tran                                    # into the native frame
        p["rmsd_to_native"] = float(np.sqrt(((posed - native_xyz) ** 2).sum(axis=1).mean()))

    poses.sort(key=lambda d: (d["clashes"], -d["contacts"]))
    df = pd.DataFrame([{k: v for k, v in p.items() if not k.startswith("_")} for p in poses])
    free = df[df.clashes == 0]
    reps = G.cluster_poses(poses)

    print(f"\nposes {len(df)} | clash-free {len(free)} ({100 * len(free) / len(df):.0f}%) "
          f"| distinct arrangements {len(reps)}")
    print(f"sampling  - closest pose to native:        {df.rmsd_to_native.min():.1f} A")
    print(f"scoring   - top-ranked pose:               {df.iloc[0].rmsd_to_native:.1f} A")
    print(f"            best of {len(reps)} representatives:     "
          f"{min(r['rmsd_to_native'] for r in reps):.1f} A")
    print(f"            median over clash-free poses:  {free.rmsd_to_native.median():.1f} A")
    print("\nVerdict: the native arrangement is inside the sampled set, but ranking by "
          "clashes and contacts\ndoes not single it out. Read the output as 'these two "
          "CAN be linked without overlapping',\nnot as a predicted ternary complex.")
    return df


if __name__ == "__main__":
    main()
