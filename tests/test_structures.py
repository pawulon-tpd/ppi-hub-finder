"""Offline tests: no network. API calls are either stubbed or built from synthetic atoms."""
import numpy as np
import pytest
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue

from ppi import geometry as G
from ppi import network
from ppi import structures as S


# ---------- ligand classification ----------

def test_ligand_kind():
    assert S.ligand_kind("EDO", 62.1) == "additive"        # cryoprotectant
    assert S.ligand_kind("ZN", 65.4) == "additive"         # ion
    assert S.ligand_kind("ATP", 507.2) == "cofactor"
    assert S.ligand_kind("XYZ", 90.0) == "additive"        # too small to be a lead
    assert S.ligand_kind("XYZ", 220.0) == "fragment"
    assert S.ligand_kind("NUT", 551.3) == "drug-like"


def test_tidy_entry_and_summary():
    raw = {"rcsb_id": "1ABC", "struct": {"title": "test"}, "exptl": [{"method": "X-RAY"}],
           "rcsb_entry_info": {"resolution_combined": [2.4, 1.8]},
           "polymer_entities": [{"rcsb_polymer_entity": {"pdbx_description": "kinase"},
                                 "rcsb_polymer_entity_container_identifiers": {
                                     "auth_asym_ids": ["A"],
                                     "reference_sequence_identifiers": [
                                         {"database_accession": "P1", "database_name": "UniProt"},
                                         {"database_accession": "X9", "database_name": "GenBank"}]}}],
           "nonpolymer_entities": [
               {"nonpolymer_comp": {"chem_comp": {"id": "LIG", "name": "lead", "formula_weight": 480}},
                "rcsb_nonpolymer_entity_container_identifiers": {"auth_asym_ids": ["A"]}},
               {"nonpolymer_comp": {"chem_comp": {"id": "GOL", "name": "glycerol",
                                                  "formula_weight": 92}},
                "rcsb_nonpolymer_entity_container_identifiers": {"auth_asym_ids": ["A"]}}]}
    e = S._tidy_entry(raw)
    assert e["resolution"] == 1.8                               # best of the list
    assert e["proteins"][0]["accessions"] == ["P1"]             # GenBank ignored
    assert [(l["id"], l["kind"]) for l in e["ligands"]] == [("LIG", "drug-like"),
                                                            ("GOL", "additive")]
    row = S.protein_summary("KIN", "P1", [e])
    assert row["n_structures"] == 1 and row["n_with_ligand"] == 1
    assert row["best_ligand_pdb"] == "1ABC" and row["n_distinct_ligands"] == 1
    assert "LIG" in row["top_ligands"]


def test_apo_only_protein_has_no_ligand_pdb():
    apo = {"pdb_id": "2XYZ", "title": "", "method": "X-RAY", "resolution": 2.0,
           "proteins": [{"accessions": ["P2"], "chains": ["A"], "description": ""}],
           "ligands": [{"id": "SO4", "name": "sulfate", "mw": 96.0, "chains": ["A"],
                        "kind": "additive"}]}
    row = S.protein_summary("APO", "P2", [apo])
    assert row["n_with_ligand"] == 0 and row["best_ligand_pdb"] == ""
    assert row["best_pdb"] == "2XYZ"


def test_cocomplex_and_pair_table():
    both = {"pdb_id": "3CO1", "title": "", "method": "X-RAY", "resolution": 2.1,
            "proteins": [{"accessions": ["PA"], "chains": ["A"], "description": ""},
                         {"accessions": ["PB"], "chains": ["B"], "description": ""}],
            "ligands": [{"id": "L1", "name": "", "mw": 400.0, "chains": ["A"],
                         "kind": "drug-like"}]}
    only_a = {"pdb_id": "4ONE", "title": "", "method": "X-RAY", "resolution": 1.5,
              "proteins": [{"accessions": ["PA"], "chains": ["A"], "description": ""}],
              "ligands": [{"id": "L2", "name": "", "mw": 350.0, "chains": ["A"],
                           "kind": "drug-like"}]}
    meta = {"3CO1": both, "4ONE": only_a}
    per = {"A": ["3CO1", "4ONE"], "B": ["3CO1"]}
    accs = {"A": "PA", "B": "PB"}
    assert [e["pdb_id"] for e in S.cocomplexes("PA", "PB", per, meta, "A", "B")] == ["3CO1"]
    row = S.pair_table([("A", "B")], accs, per, meta).iloc[0]
    assert row.n_cocomplex == 1 and row.cocomplex_pdbs == "3CO1"
    assert row.ligand_pdb_a == "4ONE"            # best resolution with a ligand
    # the co-complex ligand sits on chain A, so it must NOT count as a ligand for B
    assert not row.both_have_ligand and row.ligand_pdb_b == ""
    assert S.ligands_of(both, "PA") and not S.ligands_of(both, "PB")
    lt = S.ligand_table(meta, {"PA": "A", "PB": "B"})
    assert set(lt.ligand_id) == {"L1", "L2"} and "GOL" not in set(lt.ligand_id)


def test_accessions_skips_failures(monkeypatch):
    def fake(gene, species=9606, refresh=False):
        if gene == "BAD":
            raise RuntimeError("boom")
        return {"acc": f"P_{gene}"} if gene != "MISSING" else {}
    monkeypatch.setattr(S, "uniprot_accession", fake)
    assert S.accessions(["AAA", "BAD", "MISSING"], quiet=True) == {"AAA": "P_AAA"}


# ---------- top pairs ----------

def test_top_pairs_are_distinct_and_ordered():
    inputs = [f"I{i}" for i in range(1, 7)]
    edges = ([{"a": "HUB", "b": f"I{i}", "source": "STRING", "weight": 1} for i in (1, 2, 3, 4)]
             + [{"a": "P", "b": f"I{i}", "source": "STRING", "weight": 1} for i in (1, 2, 5)]
             + [{"a": "Q", "b": f"I{i}", "source": "STRING", "weight": 1} for i in (3, 4, 6)]
             + [{"a": "R", "b": f"I{i}", "source": "STRING", "weight": 1} for i in (5, 6)])
    G_ = network.build_graph(edges, inputs)
    table = network.reach_table(G_, inputs)
    pairs = network.top_pairs(table, inputs, n=2)
    # several pairs reach all six (HUB+R, P+Q); the best must be one of them
    assert pairs[0]["reached_pct"] == 100.0 and pairs[0]["reached_any"] == 6
    assert pairs[0]["score"] >= pairs[1]["score"]
    used = [p[k] for p in pairs for k in ("protein_a", "protein_b")]
    assert len(used) == len(set(used))                                    # distinct proteins


# ---------- geometry on synthetic atoms ----------

def _atom(name, xyz, element="C"):
    return Atom(name, np.array(xyz, "f"), 0.0, 1.0, " ", name, 1, element=element)


def _model(lig_xyz, wall_xyz):
    """A ligand plus a flat wall of protein atoms, so burial has a known direction."""
    model = Model(0)
    chain = Chain("A")
    lig = Residue((f"H_LIG", 1, " "), "LIG", " ")
    for i, p in enumerate(lig_xyz):
        lig.add(_atom(f"C{i}", p))
    chain.add(lig)
    prot = Residue((" ", 2, " "), "ALA", " ")
    for i, p in enumerate(wall_xyz):
        prot.add(_atom(f"W{i}", p))
    chain.add(prot)
    model.add(chain)
    return model, {"resname": "LIG", "chain": "A", "resseq": 1, "residue": lig}


def test_rotations():
    a, b = np.array([1.0, 0, 0]), np.array([0, 1.0, 0])
    assert np.allclose(G._rot_between(a, b) @ a, b, atol=1e-9)
    assert np.allclose(G._rot_between(a, -a) @ a, -a, atol=1e-9)          # antiparallel case
    assert np.allclose(G._rot_between(a, a), np.eye(3))
    spun = G._rot_axis([0, 0, 1.0], np.pi) @ a
    assert np.allclose(spun, -a, atol=1e-9)
    assert np.allclose(np.linalg.det(G._rot_axis([0, 1.0, 0], 0.7)), 1.0)  # proper rotation


def test_count_pairs_and_clusters():
    a = np.zeros((1, 3))
    assert G._count_pairs(a, np.array([[2.0, 0, 0]]), a[0], 0.0, 2.5, 4.5) == (1, 1)  # clash
    assert G._count_pairs(a, np.array([[4.0, 0, 0]]), a[0], 0.0, 2.5, 4.5) == (0, 1)  # contact
    assert G._count_pairs(a, np.array([[50.0, 0, 0]]), a[0], 0.0, 2.5, 4.5) == (0, 0)

    rows = [{"clashes": 0, "contacts": 50, "_centroid": np.array([0.0, 0, 0])},
            {"clashes": 0, "contacts": 40, "_centroid": np.array([1.0, 0, 0])},   # duplicate
            {"clashes": 0, "contacts": 30, "_centroid": np.array([30.0, 0, 0])},
            {"clashes": 5, "contacts": 99, "_centroid": np.array([60.0, 0, 0])}]  # clashing
    reps = G.cluster_poses(rows, cutoff=8.0)
    assert [r["contacts"] for r in reps] == [50, 30]


def test_sasa_and_exit_vector_point_away_from_the_wall():
    # ligand along x at the origin; a dense wall on the -x side
    lig = [(0.0, 0, 0), (1.5, 0, 0), (3.0, 0, 0)]
    wall = [(-4.0, y, z) for y in (-2, 0, 2) for z in (-2, 0, 2)]
    model, l = _model(lig, wall)
    s = G.ligand_sasa(model, l)
    assert 0.0 <= s["buried_fraction"] <= 1.0
    assert len(s["atoms"]) == 3
    vecs = G.exit_vectors(model, l, s, top=3)
    assert vecs and vecs[0]["direction"][0] > 0.5        # points away from the wall (+x)
    tbl = G.exit_vector_table(vecs, "TEST", "LIG")
    assert set(["pdb_id", "atom", "openness", "no_open_vector"]) <= set(tbl.columns)


def test_fully_buried_ligand_falls_back_and_is_flagged():
    lig = [(0.0, 0, 0), (1.2, 0, 0)]
    shell = [(x, y, z) for x in (-3.5, 4.7) for y in (-3.5, 0, 3.5) for z in (-3.5, 0, 3.5)]
    shell += [(0.6, y, z) for y in (-3.5, 3.5) for z in (-3.5, 0, 3.5)]
    model, l = _model(lig, shell)
    s = G.ligand_sasa(model, l)
    vecs = G.exit_vectors(model, l, s, min_rel=0.99, top=2)      # force the buried branch
    assert vecs and all(v["buried"] for v in vecs)
    assert G.exit_vector_table(vecs).no_open_vector.all()
    assert G.exit_vectors(model, l, s, min_rel=0.99, fallback=False) == []


def test_ternary_scan_geometry():
    """Two separated ligands: the scan must put B at the requested distance, facing A."""
    model_a, la = _model([(0.0, 0, 0)], [(-4.0, y, z) for y in (-2, 2) for z in (-2, 2)])
    model_b, lb = _model([(0.0, 0, 0)], [(4.0, y, z) for y in (-2, 2) for z in (-2, 2)])
    sa, sb = G.ligand_sasa(model_a, la), G.ligand_sasa(model_b, lb)
    va, vb = G.exit_vectors(model_a, la, sa, top=1), G.exit_vectors(model_b, lb, sb, top=1)
    rows = G.ternary_scan(model_a, la, va[0], model_b, lb, vb[0],
                          linkers=(10.0,), n_rot=4, tilts=(0,))
    assert len(rows) == 4
    for r in rows:
        assert r["linker_A"] == 10.0 and r["tilt_deg"] == 0.0
        assert r["clashes"] == 0                      # 10 A apart, walls facing outward
    # tilt rings add poses but keep the anchor separation
    tilted = G.ternary_scan(model_a, la, va[0], model_b, lb, vb[0],
                            linkers=(10.0,), n_rot=2, tilts=(0, 30), n_azimuth=4)
    assert len(tilted) == 2 * (1 + 4)
    assert {r["tilt_deg"] for r in tilted} == {0.0, 30.0}


def test_chains_near_ligand_ignores_packing_neighbours():
    model = Model(0)
    lig_chain = Chain("A")
    lig = Residue(("H_LIG", 1, " "), "LIG", " ")
    for i, p in enumerate([(0.0, 0, 0), (1.5, 0, 0)]):
        lig.add(_atom(f"C{i}", p))
    lig_chain.add(lig)
    real = Residue((" ", 2, " "), "ALA", " ")
    for i, p in enumerate([(3.0, y, z) for y in (-2, 0, 2) for z in (-2, 0, 2)]):
        real.add(_atom(f"R{i}", p))
    lig_chain.add(real)
    model.add(lig_chain)
    far = Chain("Z")                                  # brushes the ligand once
    res = Residue((" ", 3, " "), "ALA", " ")
    res.add(_atom("T0", (0.0, 5.5, 0.0)))
    far.add(res)
    model.add(far)
    assert G.chains_near_ligand(model, {"residue": lig}) == {"A"}


def test_pick_ligand_skips_additives():
    model, _ = _model([(0.0, 0, 0)] * 8, [(9.0, 0, 0)])
    extra = Residue(("H_GOL", 5, " "), "GOL", " ")
    for i in range(8):
        extra.add(_atom(f"G{i}", (20.0 + i, 0, 0)))
    model["A"].add(extra)
    assert G.pick_ligand(model)["resname"] == "LIG"
    with pytest.raises(ValueError):
        G.pick_ligand(model, ligand_id="NOPE")
