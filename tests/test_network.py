from ppi import network, plot


def E(a, b, src="STRING", w=0.9):
    return {"a": a, "b": b, "source": src, "weight": w}


EDGES = ([E("HUB", x) for x in "ABC"]
         + [E("A", "B", "BioGRID", 2.0), E("a", "B", w=0.95), E("D", "E", w=0.8)])
INPUTS = list("ABCDE")
CANON = {n.upper(): n for n in INPUTS}


def test_rank_and_greedy():
    G = network.build_graph(network.dedupe_edges(EDGES, CANON), INPUTS)
    assert G["A"]["B"]["sources"] == {"STRING", "BioGRID"}
    t = network.reach_table(G, INPUTS)
    r = network.rank_targets(G, INPUTS, t)
    assert r.iloc[0].protein == "HUB" and r.iloc[0].n_input_partners == 3
    picks, left = network.greedy_cover(t, INPUTS, k=3)
    assert [p["target"] for p in picks] == ["HUB", "D"] and left == []


def test_exact_beats_greedy():
    # HUB hits 1-4, P hits 1,2,5, Q hits 3,4,6: greedy takes HUB (5 total), exact finds P+Q (6)
    inputs = [f"I{i}" for i in range(1, 7)]
    edges = ([E("HUB", f"I{i}") for i in (1, 2, 3, 4)]
             + [E("P", f"I{i}") for i in (1, 2, 5)] + [E("Q", f"I{i}") for i in (3, 4, 6)])
    G = network.build_graph(edges, inputs)
    t = network.reach_table(G, inputs)
    g, _ = network.greedy_cover(t, inputs, k=2)
    x, left = network.exact_cover(t, inputs, k=2)
    assert g[-1]["reached_any"] == 5
    assert {p["target"] for p in x} == {"P", "Q"} and x[-1]["reached_any"] == 6 and not left


def test_depth2():
    # T - X - A, X - B, T - C
    inputs = ["A", "B", "C"]
    G = network.build_graph([E("T", "X"), E("X", "A"), E("X", "B"), E("T", "C")], inputs)
    t = network.reach_table(G, inputs, depth=2)
    assert t["T"] == {"C": 1, "A": 2, "B": 2}
    r = network.rank_targets(G, inputs, t, depth=2, hop_weight=0.5).set_index("protein")
    assert r.at["T", "score"] == 2.0 and r.at["X", "score"] == 2.5
    assert "A (via X)" in r.at["T", "two_hop_partners"]
    assert network.two_hop_paths(G, "T", t) == [("X", "A"), ("X", "B")]


def test_partner_counts():
    c = network.partner_counts(EDGES, INPUTS, CANON).set_index("protein")
    assert c.loc["A", "partners_any"] == 2 and c.loc["A", "partners_BioGRID"] == 1


def test_plots(tmp_path):
    inputs = ["A", "B", "C"]
    G = network.build_graph([E("T", "X"), E("X", "A"), E("X", "B"), E("T", "C")], inputs)
    t = network.reach_table(G, inputs, depth=2)
    r = network.rank_targets(G, inputs, t, depth=2)
    best = r.iloc[0].protein
    paths = network.two_hop_paths(G, best, t)
    H = plot.focus_subgraph(G, r, extra=[x for x, _ in paths] + ["T"])
    plot.draw_png(H, r, tmp_path / "n.png", best, paths)
    plot.draw_html(H, r, tmp_path / "n.html", best, paths)
    assert (tmp_path / "n.png").stat().st_size > 0 and (tmp_path / "n.html").stat().st_size > 0


def test_biogrid_synonyms(monkeypatch):
    from ppi import biogrid_api

    def rec(i, a, b, syn_b="", pub=1):
        return {"BIOGRID_INTERACTION_ID": i, "OFFICIAL_SYMBOL_A": a, "OFFICIAL_SYMBOL_B": b,
                "SYNONYMS_A": "", "SYNONYMS_B": syn_b, "EXPERIMENTAL_SYSTEM_TYPE": "physical",
                "EXPERIMENTAL_SYSTEM": "Two-hybrid", "PUBMED_ID": pub}
    raw = [rec(1, "PTEN", "PARK2", "AR-JP|PRKN"),       # old symbol -> PRKN
           rec(2, "TP53", "PRKN", pub=2),
           rec(3, "TP53", "XYZ1", "AR|PRKN")]          # ambiguous -> unchanged
    monkeypatch.setattr(biogrid_api, "fetch_raw", lambda *a, **k: raw)
    edges = biogrid_api.interactions(["PTEN", "TP53"], "key", known_names={"PTEN", "TP53", "PRKN", "AR"})
    assert {frozenset((e["a"], e["b"])) for e in edges} == {
        frozenset(("PTEN", "PRKN")), frozenset(("TP53", "PRKN")), frozenset(("TP53", "XYZ1"))}
