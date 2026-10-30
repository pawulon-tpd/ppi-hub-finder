import networkx as nx
import pandas as pd

from ppi import interactive, mechanisms as M, network


def inter(rows):
    cols = ["source", "target", "source_genesymbol", "target_genesymbol", "is_directed",
            "is_stimulation", "is_inhibition", "consensus_direction", "consensus_stimulation",
            "consensus_inhibition", "sources", "references", "type", "curation_effort"]
    out = []
    for a, b, stim, inh, *rest in rows:
        typ = rest[0] if rest else "post_translational"
        src = rest[1] if len(rest) > 1 else "SIGNOR"
        out.append([f"P_{a}", f"P_{b}", a, b, "True", str(stim), str(inh), "True", str(stim),
                    str(inh), src, "SIGNOR:1", typ, "1"])
    return pd.DataFrame(out, columns=cols)


def enz(rows):
    cols = ["enzyme", "enzyme_genesymbol", "substrate", "substrate_genesymbol", "residue_type",
            "residue_offset", "modification", "sources", "references", "curation_effort"]
    return pd.DataFrame([[f"P_{a}", a, f"P_{b}", b, res[0], res[1:], mod, src, "x:2", "1"]
                         for a, b, mod, res, src in rows], columns=cols)


def test_build_index():
    i = inter([("AKT1", "MDM2", True, False), ("MDM2", "TP53", False, True),
               ("TP53", "MDM2", True, False, "transcriptional", "CollecTRI"),
               ("AKT1", "RAF1", True, True)])
    i.loc[len(i)] = ["COMPLEX:P_A_P_B", "P_TP53", "A_B", "TP53", "True", "True", "False",
                     "True", "True", "False", "SIGNOR", "", "post_translational", "1"]
    e = enz([("AKT1", "MDM2", "phosphorylation", "S166", "PhosphoSite"),
             ("MDM2", "TP53", "polyubiquitination", "K373", "HPRD"),
             ("MDM2", "TP53", "phosphorylation", "S15", "BEL-Large-Corpus_ProtMapper;ProtMapper")])
    idx = M.build_index(i, e, ["AKT1", "MDM2", "TP53", "RAF1"])
    assert idx[("AKT1", "MDM2")]["effect"] == "activates"
    assert idx[("AKT1", "MDM2")]["kinds"] == ["phosphorylation"]
    m = idx[("MDM2", "TP53")]
    assert m["effect"] == "inhibits" and m["kinds"] == ["ubiquitination"]   # S15 text-mined: dropped
    assert not any("S15" in d for d in m["details"])
    assert idx[("TP53", "MDM2")]["kinds"] == [M.TRANSCRIPTIONAL]
    assert idx[("AKT1", "RAF1")]["effect"] == "mixed" and idx[("AKT1", "RAF1")]["kinds"] == [M.OTHER]
    assert not any("A_B" in k for pair in idx for k in pair)                # complexes skipped
    kept = M.build_index(i, e, ["MDM2", "TP53"], include_predicted=True)
    assert any("S15" in d for d in kept[("MDM2", "TP53")]["details"])


def test_chain_sign():
    assert M.chain_sign([1, -1]) == -1 and M.chain_sign([-1, -1]) == 1
    assert M.chain_sign([1, None]) is None and M.chain_sign([0, 1]) == 0
    assert M.combine_signs([1, 1, None]) == 1 and M.combine_signs([1, -1]) == 0
    assert M.combine_signs([None]) is None


def test_target_paths_prefers_informative_bridge():
    # T reaches Y via A (no mechanism) or B (T -activates-> B -inhibits-> Y); direct T-Z inhibits
    inputs = ["Y", "Z"]
    G = network.build_graph([{"a": a, "b": b, "source": "STRING", "weight": 1}
                             for a, b in [("T", "A"), ("A", "Y"), ("T", "B"), ("B", "Y"), ("T", "Z")]],
                            inputs)
    table = network.reach_table(G, inputs, depth=2)
    idx = M.build_index(inter([("T", "B", True, False), ("B", "Y", False, True),
                               ("T", "Z", False, True), ("Z", "T", True, False)]),
                        enz([]), G.nodes)
    rows = M.target_paths(G, idx, "T", table, max_bridges=1)
    by = {r["input"]: r for r in rows}
    assert by["Y"]["via"] == "B" and by["Y"]["effect_on_input"] == "inhibits"
    assert by["Z"]["effect_on_input"] == "inhibits" and by["Z"]["effect_on_target"] == "activates"
    effects, counts = M.summarize(M.target_paths(G, idx, "T", table, max_bridges=2))
    assert effects == {"Y": "inhibits", "Z": "inhibits"}       # bridge A has no data: ignored
    paths, summary, views = M.analyse_targets(G, idx, table, ["T"])
    assert summary.iloc[0].n_inhibits == 2 and set(paths.via) == {"", "A", "B"}

    data = interactive.build_data(G, network.rank_targets(G, inputs, table, 2), {"Y", "Z", "T"},
                                  views, idx, best="T")
    assert {n["id"] for n in data["nodes"]} == {"T", "A", "B", "Y", "Z"}
    assert data["views"]["T"]["pairs"] == [["A", "T"], ["A", "Y"], ["B", "T"], ["B", "Y"], ["T", "Z"]]
    assert {(m["a"], m["b"]) for m in data["mechs"]} == {("T", "B"), ("B", "Y"), ("T", "Z"), ("Z", "T")}
    assert data["kinds"] == [M.OTHER]


def test_no_mechanism_data():
    G = nx.Graph()
    inputs = ["Y"]
    G = network.build_graph([{"a": "T", "b": "Y", "source": "STRING", "weight": 1}], inputs)
    table = network.reach_table(G, inputs)
    effects, counts = M.summarize(M.target_paths(G, {}, "T", table))
    assert effects == {"Y": M.NO_DATA} and counts[M.NO_DATA] == ["Y"]


def test_direct_partner_without_mechanism_uses_bridge():
    # AKT1 binds TP53 (no directed record) but AKT1 -> MDM2 -| TP53
    inputs = ["TP53"]
    G = network.build_graph([{"a": a, "b": b, "source": "STRING", "weight": 1}
                             for a, b in [("AKT1", "TP53"), ("AKT1", "MDM2"), ("MDM2", "TP53")]],
                            inputs)
    table = network.reach_table(G, inputs, depth=2)
    idx = M.build_index(inter([("AKT1", "MDM2", True, False), ("MDM2", "TP53", False, True)]),
                        enz([]), G.nodes)
    rows = M.target_paths(G, idx, "AKT1", table)
    assert [(r["via"], r["effect_on_input"]) for r in rows] == [("", M.NO_DATA), ("MDM2", "inhibits")]
    assert M.summarize(rows)[0] == {"TP53": "inhibits"}


def test_commercial_licence_filter():
    purposes = {"SIGNOR": "commercial", "PhosphoSite": "academic",
                "HPRD": "academic", "CollecTRI": "composite"}
    assert M.commercial_ok("SIGNOR;PhosphoSite", purposes)        # one permissive source is enough
    assert M.commercial_ok("PhosphoSite_MIMP;SIGNOR_ProtMapper", purposes)
    assert not M.commercial_ok("PhosphoSite;HPRD", purposes)
    assert not M.commercial_ok("CollecTRI", purposes)             # composite cannot carry it
    assert not M.commercial_ok("NeverHeardOf", purposes)          # unlisted: cannot verify
    assert not M.commercial_ok("", purposes)

    i = inter([("A", "B", True, False, "post_translational", "SIGNOR"),
               ("B", "C", True, False, "post_translational", "PhosphoSite")])
    e = enz([("A", "B", "phosphorylation", "S1", "PhosphoSite")])
    full = M.build_index(i, e, ["A", "B", "C"])
    comm = M.build_index(i, e, ["A", "B", "C"], purposes=purposes)
    assert set(full) == {("A", "B"), ("B", "C")}
    assert set(comm) == {("A", "B")}                    # B->C had only an academic source
    assert comm[("A", "B")]["kinds"] == [M.OTHER]       # its academic phospho site is gone
