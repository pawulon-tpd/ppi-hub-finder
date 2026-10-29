"""Minimal client for the BioGRID REST service.

Docs: https://wiki.thebiogrid.org/doku.php/biogridrest
Needs a free access key from https://webservice.thebiogrid.org/
(set BIOGRID_ACCESS_KEY or pass it explicitly).
BioGRID contains only experimentally curated interactions.
"""
import os
import time
from collections import defaultdict

import requests

BASE = "https://webservice.thebiogrid.org/interactions/"
PAGE = 10000


def get_key(key=None):
    key = key or os.environ.get("BIOGRID_ACCESS_KEY")
    if not key:
        raise SystemExit("BioGRID needs an access key: register at https://webservice.thebiogrid.org/ "
                         "then set BIOGRID_ACCESS_KEY or pass --biogrid-key")
    return key


def fetch_raw(genes, access_key, tax_id=9606, include_interactors=True):
    """Raw BioGRID records involving the given gene symbols, deduplicated by interaction id.

    include_interactors=True  -> interactions with ANY partner (partner may be outside the list)
    include_interactors=False -> only interactions among the listed genes (sent as one
                                 request, since chunking would miss pairs across chunks)
    """
    records = {}
    genes = list(genes)
    size = 100 if include_interactors else max(len(genes), 1)
    for i in range(0, len(genes), size):
        chunk = genes[i:i + size]
        start = 0
        while True:
            params = {
                "accesskey": access_key, "format": "json",
                "geneList": "|".join(chunk), "searchNames": "true",
                "taxId": tax_id, "includeInteractors": str(include_interactors).lower(),
                "interSpeciesExcluded": "true", "selfInteractionsExcluded": "true",
                "start": start, "max": PAGE}
            r = requests.get(BASE, params=params, timeout=180)
            if r.status_code != 200:
                raise RuntimeError(f"BioGRID failed ({r.status_code}): {r.text[:300]}")
            data = r.json() if r.text.strip() else {}
            page = data if isinstance(data, dict) else {
                rec.get("BIOGRID_INTERACTION_ID", id(rec)): rec for rec in data}
            records.update(page)
            if len(page) < PAGE:
                break
            start += PAGE
            time.sleep(0.5)
    return list(records.values())


def _namer(known_names):
    """Map a BioGRID symbol onto a known (e.g. STRING) name.

    BioGRID sometimes keeps an older symbol (PARK2) where STRING uses the current one (PRKN).
    If the official symbol is unknown but exactly one of its synonyms is known, use that.
    """
    known = {n.upper(): n for n in known_names or ()}

    def name(symbol, synonyms):
        if not known or symbol.upper() in known:
            return known.get(symbol.upper(), symbol)
        hits = {known[s.upper()] for s in (synonyms or "").split("|") if s.upper() in known}
        return hits.pop() if len(hits) == 1 else symbol
    return name


def interactions(genes, access_key, tax_id=9606, physical_only=True,
                 min_publications=1, include_interactors=True, known_names=None):
    """Collapse BioGRID records into one edge per protein pair.

    weight = number of distinct publications supporting the pair.
    known_names: names already used elsewhere (STRING); outdated BioGRID symbols are
    mapped onto them through BioGRID's synonym lists.
    """
    name = _namer(known_names)
    pairs = defaultdict(lambda: {"systems": set(), "pubs": set(), "n": 0})
    for rec in fetch_raw(genes, access_key, tax_id, include_interactors):
        if physical_only and rec.get("EXPERIMENTAL_SYSTEM_TYPE") != "physical":
            continue
        a = name(rec["OFFICIAL_SYMBOL_A"], rec.get("SYNONYMS_A"))
        b = name(rec["OFFICIAL_SYMBOL_B"], rec.get("SYNONYMS_B"))
        if a == b:
            continue
        p = pairs[tuple(sorted((a, b)))]
        p["n"] += 1
        p["systems"].add(rec.get("EXPERIMENTAL_SYSTEM", ""))
        p["pubs"].add(rec.get("PUBMED_ID"))

    edges = []
    for (a, b), p in pairs.items():
        if len(p["pubs"]) < min_publications:
            continue
        edges.append({
            "a": a, "b": b, "source": "BioGRID", "weight": float(len(p["pubs"])),
            "n_records": p["n"], "n_publications": len(p["pubs"]),
            "systems": "; ".join(sorted(s for s in p["systems"] if s))})
    return edges
