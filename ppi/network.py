"""Graph building, reachability (depth 1 or 2) and target ranking."""
from itertools import combinations
from math import comb

import networkx as nx
import pandas as pd


# Sticky proteins that bind almost everything (chaperones, ubiquitin, 14-3-3) and
# dominate any "who binds most of my list" ranking. Used by --exclude-hubs.
PROMISCUOUS = {
    "HSP90AA1", "HSP90AB1", "HSP90B1", "HSPA1A", "HSPA1B", "HSPA4", "HSPA5", "HSPA8",
    "HSPA9", "HSPD1", "CDC37", "STUB1",
    "UBC", "UBB", "UBA52", "RPS27A",
    "YWHAB", "YWHAE", "YWHAG", "YWHAH", "YWHAQ", "YWHAZ", "SFN",
}


def drop_proteins(edges, names, keep=()):
    """Remove edges touching any protein in `names` (case-insensitive), except those in `keep`."""
    keep_up = {k.upper() for k in keep}
    drop = {n.upper() for n in names} - keep_up
    return [e for e in edges if e["a"].upper() not in drop and e["b"].upper() not in drop]


def _fixer(canon):
    canon = canon or {}
    return lambda n: canon.get(n.upper(), n)


def dedupe_edges(edges, canon=None):
    """Canonicalise names and keep one edge per (database, protein pair)."""
    fix, seen, out = _fixer(canon), set(), []
    for e in edges:
        a, b = fix(e["a"]), fix(e["b"])
        key = (e["source"], frozenset((a, b)))
        if a == b or key in seen:
            continue
        seen.add(key)
        out.append({**e, "a": a, "b": b})
    return out


def build_graph(edges, inputs, canon=None):
    """Undirected graph; edge attribute 'sources' = databases supporting the edge."""
    fix = _fixer(canon)
    G = nx.Graph()
    G.add_nodes_from(inputs)
    for e in edges:
        a, b = fix(e["a"]), fix(e["b"])
        if a == b:
            continue
        src, w = e["source"], e["weight"]
        if G.has_edge(a, b):
            G[a][b]["sources"].add(src)
            G[a][b][f"w_{src}"] = max(G[a][b].get(f"w_{src}", 0), w)
        else:
            G.add_edge(a, b, sources={src}, **{f"w_{src}": w})
    inp = set(inputs)
    for n in G.nodes:
        G.nodes[n]["is_input"] = n in inp
    return G


def partner_counts(edges, inputs, canon=None):
    """How many distinct partners each input protein has, per database."""
    fix = _fixer(canon)
    inp = set(inputs)
    per = {n: {} for n in inputs}
    for e in edges:
        a, b = fix(e["a"]), fix(e["b"])
        for x, y in ((a, b), (b, a)):
            if x in inp and x != y:
                per[x].setdefault(e["source"], set()).add(y)
    sources = sorted({e["source"] for e in edges})
    rows = []
    for n in inputs:
        row, allp = {"protein": n}, set()
        for s in sources:
            ps = per[n].get(s, set())
            row[f"partners_{s}"] = len(ps)
            allp |= ps
        row["partners_any"] = len(allp)
        rows.append(row)
    return pd.DataFrame(rows).sort_values("partners_any", ascending=False)


# ---------- reachability ----------

def reach_table(G, inputs, depth=1):
    """{protein: {input protein: hops}} for every protein in the graph.

    hops 0 = the protein itself is an input, 1 = direct partner,
    2 = reached through one intermediate protein (only when depth == 2).
    """
    inp = set(inputs)
    nb = {x: set(G.neighbors(x)) & inp for x in G.nodes}
    table = {}
    for n in G.nodes:
        r = {}
        if depth >= 2:
            for x in G.neighbors(n):
                for y in nb[x]:
                    if y != n:
                        r[y] = 2
        for y in nb[n]:
            r[y] = 1
        if n in inp:
            r[n] = 0
        table[n] = r
    return table


def _w(hops, hop_weight):
    return 1.0 if hops <= 1 else hop_weight


def _value(r, hop_weight):
    return sum(_w(h, hop_weight) for h in r.values())


def bridges(G, n, y, max_n=3):
    """Intermediate proteins on n -> X -> y paths."""
    return sorted(x for x in G.neighbors(n) if G.has_edge(x, y))[:max_n]


def two_hop_paths(G, n, table, max_paths=40):
    """(bridge, input) pairs for inputs that n reaches only in 2 hops (for drawing)."""
    out = []
    for y, h in sorted(table[n].items()):
        if h == 2:
            out += [(x, y) for x in bridges(G, n, y, max_n=2)]
    return out[:max_paths]


def rank_targets(G, inputs, table=None, depth=1, hop_weight=0.5, annotate_top=50):
    """Rank every protein by how much of the input list it reaches.

    score = direct input partners + hop_weight * inputs reached in 2 hops.
    """
    table = table or reach_table(G, inputs, depth)
    rows = []
    for n, r in table.items():
        direct = sorted(i for i, h in r.items() if h == 1)
        two = sorted(i for i, h in r.items() if h == 2)
        rows.append({
            "protein": n, "in_input": G.nodes[n]["is_input"],
            "n_input_partners": len(direct), "n_two_hop": len(two),
            "score": round(len(direct) + hop_weight * len(two), 2),
            "coverage": len(r), "coverage_pct": round(100 * len(r) / len(inputs), 1),
            "input_partners": ", ".join(direct), "two_hop_partners": ""})
    df = (pd.DataFrame(rows)
          .sort_values(["score", "coverage", "protein"], ascending=[False, False, True])
          .reset_index(drop=True))
    if depth >= 2:
        for idx in df.index[:annotate_top]:
            n = df.at[idx, "protein"]
            two = [y for y, h in table[n].items() if h == 2]
            df.at[idx, "two_hop_partners"] = ", ".join(
                f"{y} (via {'/'.join(bridges(G, n, y))})" for y in sorted(two))
    return df


# ---------- choosing several targets ----------

def _summarize(combo, table, inputs, hop_weight):
    inp, cur, best_h, picks = set(inputs), {}, {}, []
    for t in combo:
        newly = []
        for i, h in sorted(table[t].items()):
            w = _w(h, hop_weight)
            if w > cur.get(i, 0.0):
                newly.append(i if h <= 1 else f"{i}(2-hop)")
                cur[i] = w
            best_h[i] = min(best_h.get(i, 9), h)
        picks.append({
            "target": t, "in_input": t in inp,
            "newly_reached": ", ".join(newly),
            "reached_direct": sum(1 for h in best_h.values() if h <= 1),
            "reached_any": len(cur),
            "reached_pct": round(100 * len(cur) / len(inp), 1),
            "total_score": round(sum(cur.values()), 2)})
    return picks, sorted(inp - set(cur))


def greedy_cover(table, inputs, k=3, hop_weight=0.5):
    """Pick k DIFFERENT targets one at a time, each adding the most not-yet-reached inputs.
    An input reached by several targets counts once (at its best hop distance)."""
    cur, chosen = {}, []
    nodes = sorted(table)
    for _ in range(k):
        best, gain = None, 1e-12
        for n in nodes:
            if n in chosen:
                continue
            g = sum(max(0.0, _w(h, hop_weight) - cur.get(i, 0.0)) for i, h in table[n].items())
            if g > gain:
                best, gain = n, g
        if best is None:
            break
        chosen.append(best)
        for i, h in table[best].items():
            cur[i] = max(cur.get(i, 0.0), _w(h, hop_weight))
    return _summarize(chosen, table, inputs, hop_weight)


def top_pairs(table, inputs, n=3, hop_weight=0.5, candidates=40, distinct=True):
    """The n best pairs of targets by combined weighted coverage of the input list.

    `--cover 2` returns only the single best pair; for structural follow-up it helps to
    see several. With distinct=True a protein appears in at most one returned pair, so
    the list is not n variations on the same dominant hub.
    """
    cand = sorted(table, key=lambda x: (-_value(table[x], hop_weight), x))[:candidates]
    scored = []
    for a, b in combinations(cand, 2):
        cur = {}
        for t in (a, b):
            for i, h in table[t].items():
                w = _w(h, hop_weight)
                if w > cur.get(i, 0.0):
                    cur[i] = w
        scored.append((round(sum(cur.values()), 2), len(cur), a, b))
    scored.sort(key=lambda t: (-t[0], -t[1], t[2], t[3]))
    out, used = [], set()
    for value, reached, a, b in scored:
        if distinct and (a in used or b in used):
            continue
        out.append({"protein_a": a, "protein_b": b, "score": value,
                    "reached_any": reached,
                    "reached_pct": round(100 * reached / len(inputs), 1)})
        used.update((a, b))
        if len(out) >= n:
            break
    return out


def exact_cover(table, inputs, k=3, hop_weight=0.5, candidates=60):
    """Test every combination of k targets among the `candidates` best single targets."""
    cand = sorted(table, key=lambda n: (-_value(table[n], hop_weight), n))[:candidates]
    k = min(k, len(cand))
    if comb(len(cand), k) > 5_000_000:
        raise SystemExit(f"{comb(len(cand), k):,} combinations - lower --candidates or --cover")
    best, best_v = (), -1.0
    for combo in combinations(cand, k):
        cur = {}
        for t in combo:
            for i, h in table[t].items():
                w = _w(h, hop_weight)
                if w > cur.get(i, 0.0):
                    cur[i] = w
        v = sum(cur.values())
        if v > best_v + 1e-12:
            best, best_v = combo, v
    best = sorted(best, key=lambda n: (-_value(table[n], hop_weight), n))
    return _summarize(best, table, inputs, hop_weight)
