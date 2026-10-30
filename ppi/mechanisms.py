"""How proteins act on each other: direction, sign and mechanism from OmniPath.

OmniPath (https://omnipathdb.org) combines ~100 resources (SIGNOR, PhosphoSitePlus,
KEGG, Reactome, CollecTRI, ...) into
  - interactions: directed, signed (activation / inhibition), post-translational or
    transcriptional,
  - enz_sub: enzyme -> substrate modifications with residues (phosphorylation S166, ...).

Both full tables (~40 MB for human) are downloaded once and cached. All lookups are local.
"""
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
import requests

BASE = "https://omnipathdb.org"
CACHE = Path.home() / ".cache" / "ppi-hub-finder"
MAX_AGE_DAYS = 30
SPECIES = {9606, 10090, 10116}

# Records supported ONLY by these (text mining, kinase-substrate predictions) are dropped
# unless include_predicted=True. "ProtMapper" alone is an aggregation layer, not evidence.
LOW_CONFIDENCE = ("REACH", "Sparser", "RLIMS-P", "MIMP", "NetworKIN", "BEL-Large-Corpus")

# Normalise modification names to a few toggleable kinds
KIND_ALIASES = {
    "polyubiquitination": "ubiquitination", "monoubiquitination": "ubiquitination",
    "trimethylation": "methylation", "post translational modification": "other modification",
}
TRANSCRIPTIONAL = "transcriptional regulation"
OTHER = "other regulation"

EFFECT_LABEL = {1: "activates", -1: "inhibits", 0: "mixed", None: "sign unknown"}
NO_DATA = "no mechanism data"


# ---------- download / cache ----------

def _download(endpoint, params, path):
    r = requests.get(f"{BASE}/{endpoint}", params=params, timeout=600)
    if r.status_code != 200:
        raise RuntimeError(f"OmniPath {endpoint} failed ({r.status_code}): {r.text[:300]}")
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".part")
    tmp.write_bytes(r.content)
    tmp.replace(path)


def load_tables(species=9606, license="academic", refresh=False, cache_dir=CACHE):
    """(interactions, enz_sub) DataFrames for one organism, cached for MAX_AGE_DAYS."""
    if species not in SPECIES:
        raise ValueError(f"OmniPath covers taxon ids {sorted(SPECIES)}, not {species}")
    cache_dir = Path(cache_dir)
    jobs = {
        "interactions": {"datasets": "omnipath,collectri",
                         "fields": "sources,references,type,curation_effort"},
        "enz_sub": {"fields": "sources,references,curation_effort"},
    }
    out = []
    for endpoint, extra in jobs.items():
        path = cache_dir / f"omnipath_{endpoint}_{species}_{license}.tsv"
        stale = not path.exists() or time.time() - path.stat().st_mtime > MAX_AGE_DAYS * 86400
        if refresh or stale:
            print(f"  downloading OmniPath {endpoint} (one-off, cached in {cache_dir})...")
            _download(endpoint, {"genesymbols": "yes", "organisms": species,
                                 "license": license, **extra}, path)
        out.append(pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False))
    return tuple(out)


# ---------- building the directed mechanism index ----------

def _confident(sources, include_predicted):
    if include_predicted:
        return True
    return any(s != "ProtMapper" and not any(t in s for t in LOW_CONFIDENCE)
               for s in sources.split(";") if s)


def _sign(stim, inhib):
    if stim and inhib:
        return 0
    if stim:
        return 1
    if inhib:
        return -1
    return None


def _row_sign(r):
    t = lambda k: r.get(k) == "True"  # noqa: E731
    s = _sign(t("consensus_stimulation"), t("consensus_inhibition"))
    return s if s is not None else _sign(t("is_stimulation"), t("is_inhibition"))


def combine_signs(signs):
    """Several records for the same step: agree -> that sign, disagree -> mixed (0)."""
    known = {s for s in signs if s is not None}
    if not known:
        return None
    return known.pop() if len(known) == 1 else 0


def build_index(interactions, enz_sub, names, include_predicted=False):
    """{(source, target): mechanism} for ordered pairs of proteins in `names`.

    mechanism = {"sign": 1 / -1 / 0 (mixed) / None (directed, sign unknown),
                 "effect": label, "kinds": sorted kinds, "details": [lines],
                 "sources": "...", "n_refs": int}
    Names are matched case-insensitively and returned in the spelling used in `names`.
    """
    canon = {n.upper(): n for n in names}
    acc = defaultdict(lambda: {"signs": [], "kinds": set(), "details": [], "sources": set(),
                               "refs": set(), "sites": defaultdict(list)})

    for r in interactions.to_dict("records"):
        if r["source"].startswith("COMPLEX:") or r["target"].startswith("COMPLEX:"):
            continue
        a = canon.get(r["source_genesymbol"].upper())
        b = canon.get(r["target_genesymbol"].upper())
        if not a or not b or a == b or r.get("is_directed") != "True":
            continue
        if not _confident(r["sources"], include_predicted):
            continue
        sign = _row_sign(r)
        kind = TRANSCRIPTIONAL if r["type"] == "transcriptional" else None
        m = acc[(a, b)]
        m["signs"].append(sign)
        if kind:
            m["kinds"].add(kind)
            m["details"].append(f"transcriptional: {EFFECT_LABEL[sign]}")
        else:
            m["details"].append(f"post-translational: {EFFECT_LABEL[sign]}")
            m["post"] = True
        m["sources"].update(s.split("_")[0] for s in r["sources"].split(";") if s)
        m["refs"].update(x.split(":")[-1] for x in r["references"].split(";") if x)

    for r in enz_sub.to_dict("records"):
        a = canon.get(r["enzyme_genesymbol"].upper())
        b = canon.get(r["substrate_genesymbol"].upper())
        if not a or not b or a == b or not _confident(r["sources"], include_predicted):
            continue
        kind = KIND_ALIASES.get(r["modification"], r["modification"])
        m = acc[(a, b)]
        m["kinds"].add(kind)
        m["sites"][kind].append(f"{r['residue_type']}{r['residue_offset']}")
        m["sources"].update(s.split("_")[0] for s in r["sources"].split(";") if s)
        m["refs"].update(x.split(":")[-1] for x in r["references"].split(";") if x)

    index = {}
    for pair, m in acc.items():
        if m.get("post") and not m["sites"]:
            m["kinds"].add(OTHER)
        details = sorted(set(m["details"]))
        for kind, sites in sorted(m["sites"].items()):
            uniq = sorted(set(sites), key=lambda s: (len(s), s))
            details.append(f"{kind} {', '.join(uniq[:8])}{' ...' if len(uniq) > 8 else ''}")
        sign = combine_signs(m["signs"])
        index[pair] = {"sign": sign, "effect": EFFECT_LABEL[sign], "kinds": sorted(m["kinds"]),
                       "details": details, "sources": ", ".join(sorted(m["sources"])),
                       "n_refs": len(m["refs"])}
    return index


def index_table(index, G=None):
    """DataFrame of all directed mechanisms (for mechanisms.csv)."""
    rows = [{"source": a, "target": b, "effect": m["effect"], "kinds": "; ".join(m["kinds"]),
             "details": "; ".join(m["details"]), "n_references": m["n_refs"],
             "databases": m["sources"],
             **({"in_ppi_network": G.has_edge(a, b)} if G is not None else {})}
            for (a, b), m in index.items()]
    cols = ["source", "target", "effect", "kinds", "details", "n_references", "databases"]
    return pd.DataFrame(rows, columns=cols + (["in_ppi_network"] if G is not None else []))


# ---------- paths and per-target effects ----------

def chain_sign(steps):
    """Sign of a chain of directed steps: product of signs; any unknown -> None, any mixed -> 0."""
    if any(s is None for s in steps):
        return None
    out = 1
    for s in steps:
        out *= s
    return out


def _chain(index, nodes):
    """Mechanisms along consecutive node pairs; None if any step has no directed record."""
    steps = [index.get((x, y)) for x, y in zip(nodes, nodes[1:])]
    return steps if all(steps) else None


def _describe(index, nodes):
    steps = _chain(index, nodes)
    if steps is None:
        return None, ""
    sign = chain_sign([m["sign"] for m in steps])
    text = " | ".join(f"{x}->{y}: {m['effect']} ({'; '.join(m['details'])})"
                      for (x, y), m in zip(zip(nodes, nodes[1:]), steps))
    return sign, text


def _bridge_rank(index, t, x, y):
    """Sort key for bridges: signed forward chain < any forward chain < reverse < nothing;
    then clear sign before mixed; then best-supported chain (papers behind its weakest step)."""
    f = _chain(index, [t, x, y])
    r = _chain(index, [y, x, t])
    fs = chain_sign([m["sign"] for m in f]) if f else None
    rs = chain_sign([m["sign"] for m in r]) if r else None
    cat = 0 if f and fs is not None else 1 if f else 2 if r and rs is not None else 3 if r else 4
    best = f or r
    sign = fs if f else rs
    support = min(m["n_refs"] for m in best) if best else 0
    return (cat, sign == 0, -support, x)


def target_paths(G, index, target, table, max_bridges=3):
    """One row per (input reached, bridge): how `target` acts on the input and vice versa.

    forward = target -> (bridge ->) input; reverse = input -> (bridge ->) target.
    A direct partner with no directed record also gets bridges that DO give a directed
    chain (e.g. AKT1 binds TP53, but acts on it through AKT1 -> MDM2 -| TP53).
    """
    rows = []
    for y, hops in sorted(table[target].items()):
        if hops == 0:
            continue
        cand = sorted((x for x in G.neighbors(target) if G.has_edge(x, y)),
                      key=lambda x: _bridge_rank(index, target, x, y))
        if hops == 1:
            vias = [None]
            if (target, y) not in index:
                vias += [x for x in cand if _bridge_rank(index, target, x, y)[0] <= 1][:max_bridges]
        else:
            vias = cand[:max_bridges]
        for x in vias:
            nodes = [target, y] if x is None else [target, x, y]
            fs, ft = _describe(index, nodes)
            rs, rt = _describe(index, nodes[::-1])
            rows.append({
                "target": target, "input": y, "hops": hops, "via": x or "",
                "path": " - ".join(nodes),
                "effect_on_input": EFFECT_LABEL[fs] if ft else NO_DATA, "how": ft,
                "effect_on_target": EFFECT_LABEL[rs] if rt else NO_DATA, "how_reverse": rt,
                "_sign": fs, "_has": bool(ft)})
    return rows


def summarize(rows):
    """Per input: overall effect of the target on it (bridges combined), then counts."""
    per = defaultdict(list)
    for r in rows:
        per[r["input"]].append(r)
    effects = {}
    for y, rs in per.items():
        known = [r["_sign"] for r in rs if r["_has"]]
        if not known:
            effects[y] = NO_DATA
        else:
            effects[y] = EFFECT_LABEL[combine_signs(known)] if any(
                s is not None for s in known) else EFFECT_LABEL[None]
    counts = {lab: sorted(y for y, e in effects.items() if e == lab)
              for lab in ("activates", "inhibits", "mixed", "sign unknown", NO_DATA)}
    return effects, counts


def analyse_targets(G, index, table, targets, max_bridges=3):
    """(paths DataFrame, per-target summary DataFrame, {target: rows}) for the given targets."""
    all_rows, summary, views = [], [], {}
    for t in targets:
        rows = target_paths(G, index, t, table, max_bridges)
        views[t] = rows
        all_rows += rows
        _, counts = summarize(rows)
        summary.append({"protein": t, **{f"n_{k.replace(' ', '_')}": len(v)
                                         for k, v in counts.items()},
                        **{k.replace(" ", "_"): ", ".join(v) for k, v in counts.items()
                           if k != NO_DATA}})
    paths = pd.DataFrame(all_rows).drop(columns=["_sign", "_has"], errors="ignore")
    return paths, pd.DataFrame(summary), views
