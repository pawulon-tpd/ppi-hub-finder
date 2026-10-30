"""What is known structurally about a protein: PDB entries, bound ligands, co-complexes.

Sources (both free, no key):
  - UniProt REST   gene symbol -> reviewed accession
  - RCSB search + data (GraphQL) APIs: entries per accession, resolution, method,
    chains, and every bound small molecule.

Everything is cached under ~/.cache/ppi-hub-finder/rcsb/ so reruns cost no API calls.
"""
import json
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
import requests

from .mechanisms import CACHE

UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL = "https://data.rcsb.org/graphql"
PAUSE = 0.2
BATCH = 50
MAX_AGE_DAYS = 30

# Crystallisation additives, buffers, cryoprotectants and ions: present in the crystal,
# not a chemical starting point. Anything here is reported but never counted as a ligand.
ADDITIVES = {
    "HOH", "DOD", "EDO", "GOL", "PEG", "PGE", "PG4", "PG0", "1PE", "2PE", "P6G", "7PE",
    "MPD", "BU3", "IPA", "MOH", "ETH", "DMS", "DMF", "ACN", "ACT", "ACY", "FMT", "CIT",
    "FLC", "TLA", "MLA", "MES", "EPE", "TRS", "BTB", "BIS", "CAC", "IMD", "TAR", "SIN",
    "SO4", "PO4", "NO3", "CO3", "BCT", "CL", "BR", "IOD", "F", "NA", "K", "LI", "CS",
    "RB", "MG", "CA", "ZN", "MN", "FE", "FE2", "CU", "CU1", "NI", "CO", "CD", "HG",
    "PT", "AU", "AG", "SR", "BA", "PB", "AL", "GA", "SE", "MSE", "UNX", "UNL", "UNK",
    "BME", "DTT", "DTU", "TCE", "MRD", "SCN", "AZI", "NH4", "OXY", "PER", "PEO",
    "NAG", "NDG", "BMA", "MAN", "GAL", "GLC", "FUC", "XYS", "SIA", "BGC", "A2G",
    "PLM", "MYR", "OLA", "STE", "LDA", "LMT", "DDQ", "C8E", "BOG", "OCT", "HEX",
}

# Biologically meaningful but not a drug-like starting point on their own. For a kinase the
# ATP site IS the pocket you want, so these are counted separately rather than discarded.
COFACTORS = {
    "ATP", "ADP", "AMP", "ANP", "ACP", "AGS", "APC", "ADX", "GTP", "GDP", "GMP", "GNP",
    "GSP", "CTP", "UTP", "TTP", "NAD", "NAI", "NAP", "NDP", "FAD", "FMN", "FDA",
    "SAM", "SAH", "COA", "ACO", "HEM", "HEC", "HEA", "BCL", "CLA", "PLP", "TPP",
    "BTN", "B12", "F43", "MQ7", "UQ1", "RET", "ATR",
}

DRUGLIKE_MIN_MW = 150.0   # below this a "ligand" is an ion or a fragment of solvent
FRAGMENT_MAX_MW = 300.0   # 150-300 = fragment-sized; above = lead/drug-sized


RETRIES = 4


def _request(method, url, **kw):
    """One HTTP call with retries on transient failures (5xx, 429, timeouts)."""
    last = None
    for attempt in range(RETRIES):
        try:
            r = requests.request(method, url, timeout=kw.pop("timeout", 120), **kw)
            if r.status_code < 400:
                time.sleep(PAUSE)
                return r
            if r.status_code not in (429, 500, 502, 503, 504):
                r.raise_for_status()
            last = RuntimeError(f"{url} -> {r.status_code}")
        except requests.RequestException as exc:
            last = exc
        time.sleep(2 ** attempt * PAUSE)
    raise last


# ---------- tiny cached-JSON helper ----------

def _cache_path(*parts):
    p = Path(CACHE) / "rcsb"
    p.mkdir(parents=True, exist_ok=True)
    return p.joinpath(*parts)


def _cached(name, build, refresh=False, max_age_days=MAX_AGE_DAYS):
    path = _cache_path(name)
    if not refresh and path.exists():
        if time.time() - path.stat().st_mtime < max_age_days * 86400:
            return json.loads(path.read_text())
    data = build()
    path.write_text(json.dumps(data))
    return data


# ---------- identifiers ----------

def uniprot_accession(gene, species=9606, refresh=False):
    """Reviewed (Swiss-Prot) accession for a gene symbol, or None."""
    def build():
        r = _request("GET", UNIPROT, timeout=60, params={
            "query": f"gene_exact:{gene} AND organism_id:{species} AND reviewed:true",
            "fields": "accession,gene_primary,protein_name", "format": "tsv", "size": 5})
        rows = [ln.split("\t") for ln in r.text.strip().splitlines()[1:]]
        return {"acc": rows[0][0], "name": rows[0][2]} if rows else {}
    return _cached(f"uniprot_{species}_{gene.upper()}.json", build, refresh)


def accessions(genes, species=9606, refresh=False, quiet=False):
    """{gene: accession} for the genes that have a reviewed entry."""
    out = {}
    for g in genes:
        try:
            hit = uniprot_accession(g, species, refresh)
        except Exception as exc:                                  # noqa: BLE001
            if not quiet:
                print(f"  UniProt lookup failed for {g}: {exc}")
            continue
        if hit.get("acc"):
            out[g] = hit["acc"]
    return out


# ---------- RCSB ----------

def entries_for_accession(acc, refresh=False):
    """PDB entry IDs whose polymer entities reference this UniProt accession."""
    def build():
        q = {"query": {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers."
                             "reference_sequence_identifiers.database_accession",
                "operator": "exact_match", "value": acc}},
             "return_type": "entry",
             "request_options": {"paginate": {"start": 0, "rows": 2000},
                                 "results_content_type": ["experimental"]}}
        r = _request("POST", SEARCH, json=q)
        if r.status_code == 204:                                   # no hits
            return []
        return [h["identifier"] for h in r.json().get("result_set", [])]
    return _cached(f"entries_{acc}.json", build, refresh)


ENTRY_QUERY = """query($ids:[String!]!){entries(entry_ids:$ids){
  rcsb_id struct{title} exptl{method}
  rcsb_entry_info{resolution_combined}
  polymer_entities{
    rcsb_polymer_entity{pdbx_description}
    rcsb_polymer_entity_container_identifiers{
      auth_asym_ids
      reference_sequence_identifiers{database_accession database_name}}}
  nonpolymer_entities{
    nonpolymer_comp{chem_comp{id name formula_weight}}
    rcsb_nonpolymer_entity_container_identifiers{auth_asym_ids}}}}"""


def _fetch_entries(ids):
    r = _request("POST", GRAPHQL, timeout=180,
                 json={"query": ENTRY_QUERY, "variables": {"ids": list(ids)}})
    body = r.json()
    if body.get("errors"):
        raise RuntimeError(f"RCSB GraphQL: {body['errors'][0].get('message')}")
    return [e for e in body.get("data", {}).get("entries") or [] if e]


def entry_metadata(entry_ids, refresh=False, quiet=False):
    """{entry_id: tidy metadata dict}, fetched in batches and cached one file per entry."""
    out, todo = {}, []
    for eid in dict.fromkeys(entry_ids):
        path = _cache_path(f"entry_{eid}.json")
        if not refresh and path.exists():
            out[eid] = json.loads(path.read_text())
        else:
            todo.append(eid)
    if todo and not quiet:
        print(f"  fetching {len(todo)} PDB entries from RCSB...")
    for i in range(0, len(todo), BATCH):
        for raw in _fetch_entries(todo[i:i + BATCH]):
            tidy = _tidy_entry(raw)
            _cache_path(f"entry_{tidy['pdb_id']}.json").write_text(json.dumps(tidy))
            out[tidy["pdb_id"]] = tidy
    return out


def _tidy_entry(e):
    res = (e.get("rcsb_entry_info") or {}).get("resolution_combined") or []
    prot = []
    for pe in e.get("polymer_entities") or []:
        ids = pe.get("rcsb_polymer_entity_container_identifiers") or {}
        accs = [r["database_accession"] for r in (ids.get("reference_sequence_identifiers") or [])
                if r.get("database_name") == "UniProt"]
        prot.append({"accessions": accs, "chains": ids.get("auth_asym_ids") or [],
                     "description": (pe.get("rcsb_polymer_entity") or {}).get("pdbx_description")})
    lig = []
    for ne in e.get("nonpolymer_entities") or []:
        cc = ((ne.get("nonpolymer_comp") or {}).get("chem_comp")) or {}
        if not cc.get("id"):
            continue
        lig.append({"id": cc["id"], "name": cc.get("name") or "",
                    "mw": float(cc.get("formula_weight") or 0.0),
                    "chains": (ne.get("rcsb_nonpolymer_entity_container_identifiers")
                               or {}).get("auth_asym_ids") or [],
                    "kind": ligand_kind(cc["id"], float(cc.get("formula_weight") or 0.0))})
    return {"pdb_id": e["rcsb_id"], "title": (e.get("struct") or {}).get("title") or "",
            "method": ((e.get("exptl") or [{}])[0] or {}).get("method") or "",
            "resolution": min(res) if res else None, "proteins": prot, "ligands": lig}


def ligand_kind(comp_id, mw):
    """'drug-like' | 'fragment' | 'cofactor' | 'additive'."""
    if comp_id in ADDITIVES:
        return "additive"
    if comp_id in COFACTORS:
        return "cofactor"
    if mw < DRUGLIKE_MIN_MW:
        return "additive"
    return "fragment" if mw < FRAGMENT_MAX_MW else "drug-like"


# ---------- per-protein summary ----------

REAL_LIGAND = ("drug-like", "fragment")


def chains_of(entry, acc):
    """auth chain ids in this entry that belong to the given UniProt accession."""
    return {c for p in entry["proteins"] if acc in p["accessions"] for c in p["chains"]}


def ligands_of(entry, acc=None):
    """Real ligands in this entry; with `acc`, only those on that protein's own chains.

    In a co-complex the entry has a ligand, but it may sit on the partner. Nonpolymer
    entities carry the auth chain of the protein they bind, so the chain sets decide.
    """
    out = [l for l in entry["ligands"] if l["kind"] in REAL_LIGAND]
    if acc is None:
        return out
    mine = chains_of(entry, acc)
    return [l for l in out if not mine or not l["chains"] or (set(l["chains"]) & mine)]


def _best(entries, want_ligand, acc=None):
    """Highest-resolution entry, optionally restricted to ones with a real ligand."""
    pool = [e for e in entries if not want_ligand or ligands_of(e, acc)]
    pool = [e for e in pool if e["resolution"] is not None] or pool
    return min(pool, key=lambda e: (e["resolution"] is None, e["resolution"] or 9e9,
                                    e["pdb_id"]), default=None)


def protein_summary(gene, acc, entries):
    """One row: how much structural information exists for this protein."""
    ligs = {}
    for e in entries:
        for l in ligands_of(e, acc):
            prev = ligs.get(l["id"])
            if prev is None or (e["resolution"] or 9e9) < prev[1]:
                ligs[l["id"]] = (l, e["resolution"] or 9e9, e["pdb_id"])
    with_lig = [e for e in entries if ligands_of(e, acc)]
    best_any, best_lig = _best(entries, False), _best(entries, True, acc)
    ranked = sorted(ligs.values(), key=lambda t: (-t[0]["mw"], t[2]))
    return {
        "protein": gene, "uniprot": acc, "n_structures": len(entries),
        "n_with_ligand": len(with_lig),
        "best_resolution": best_any["resolution"] if best_any else None,
        "best_pdb": best_any["pdb_id"] if best_any else "",
        "best_ligand_pdb": best_lig["pdb_id"] if best_lig else "",
        "best_ligand_resolution": best_lig["resolution"] if best_lig else None,
        "n_distinct_ligands": len(ligs),
        "top_ligands": ", ".join(f"{l['id']} ({l['mw']:.0f} Da, {pdb})"
                                 for l, _, pdb in ranked[:5]),
        "example_pdbs": ", ".join(e["pdb_id"] for e in sorted(
            with_lig or entries,
            key=lambda e: (e["resolution"] is None, e["resolution"] or 9e9))[:5]),
    }


def collect(genes, species=9606, refresh=False, quiet=False):
    """(summary DataFrame, {gene: accession}, {entry_id: metadata}, {gene: [entry_ids]})."""
    accs = accessions(genes, species, refresh, quiet)
    per_gene, all_ids = {}, []
    for g, acc in accs.items():
        ids = entries_for_accession(acc, refresh)
        per_gene[g] = ids
        all_ids += ids
    meta = entry_metadata(all_ids, refresh, quiet)
    rows = [protein_summary(g, accs[g], [meta[i] for i in ids if i in meta])
            for g, ids in per_gene.items()]
    df = (pd.DataFrame(rows).sort_values(["n_with_ligand", "n_structures"], ascending=False)
          if rows else pd.DataFrame())
    return df, accs, meta, per_gene


# ---------- pairs ----------

def cocomplexes(acc_a, acc_b, per_gene_ids, meta, gene_a, gene_b):
    """Entries containing BOTH proteins, best resolution first."""
    shared = set(per_gene_ids.get(gene_a, [])) & set(per_gene_ids.get(gene_b, []))
    out = []
    for eid in shared:
        e = meta.get(eid)
        if not e:
            continue
        accs_in = {a for p in e["proteins"] for a in p["accessions"]}
        if acc_a in accs_in and acc_b in accs_in:
            out.append(e)
    return sorted(out, key=lambda e: (e["resolution"] is None, e["resolution"] or 9e9))


def pair_table(pairs, accs, per_gene_ids, meta):
    """One row per (target pair): structures each side, and any structure holding both."""
    rows = []
    for a, b in pairs:
        ea = [meta[i] for i in per_gene_ids.get(a, []) if i in meta]
        eb = [meta[i] for i in per_gene_ids.get(b, []) if i in meta]
        co = (cocomplexes(accs[a], accs[b], per_gene_ids, meta, a, b)
              if a in accs and b in accs else [])
        acc_a, acc_b = accs.get(a), accs.get(b)
        la, lb = _best(ea, True, acc_a), _best(eb, True, acc_b)
        rows.append({
            "protein_a": a, "protein_b": b,
            "ligand_pdb_a": la["pdb_id"] if la else "", "ligand_pdb_b": lb["pdb_id"] if lb else "",
            "resolution_a": la["resolution"] if la else None,
            "resolution_b": lb["resolution"] if lb else None,
            "ligands_a": ", ".join(sorted({l["id"] for l in ligands_of(la, acc_a)})) if la else "",
            "ligands_b": ", ".join(sorted({l["id"] for l in ligands_of(lb, acc_b)})) if lb else "",
            "n_cocomplex": len(co),
            "cocomplex_pdbs": ", ".join(e["pdb_id"] for e in co[:5]),
            "best_cocomplex_resolution": co[0]["resolution"] if co else None,
            "both_have_ligand": bool(la and lb)})
    return pd.DataFrame(rows)


def ligand_table(meta, genes_by_acc):
    """Every non-additive ligand seen, with the protein and entry it came from."""
    rows = []
    for e in meta.values():
        accs_in = {a for p in e["proteins"] for a in p["accessions"]}
        who = sorted({genes_by_acc[a] for a in accs_in if a in genes_by_acc})
        for l in e["ligands"]:
            if l["kind"] == "additive":
                continue
            rows.append({"pdb_id": e["pdb_id"], "proteins": ", ".join(who),
                         "ligand_id": l["id"], "kind": l["kind"], "mw": round(l["mw"], 1),
                         "chains": ", ".join(l["chains"]), "ligand_name": l["name"][:120],
                         "resolution": e["resolution"], "method": e["method"],
                         "title": e["title"][:120]})
    df = pd.DataFrame(rows)
    return df.sort_values(["proteins", "kind", "mw"], ascending=[True, True, False]) \
        if len(df) else df
