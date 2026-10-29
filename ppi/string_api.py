"""Minimal client for the STRING REST API (https://string-db.org/help/api/).

No key needed. STRING asks callers to wait ~1 s between requests and to send
a caller_identity, both handled here.
"""
import time

import requests

BASE = "https://string-db.org/api/json"
CALLER = "ppi-hub-finder"
PAUSE = 1.0
EVIDENCE = ["nscore", "fscore", "pscore", "ascore", "escore", "dscore", "tscore"]


def _post(endpoint, params, retries=4):
    params = {**params, "caller_identity": CALLER}
    for attempt in range(retries):
        r = requests.post(f"{BASE}/{endpoint}", data=params, timeout=120)
        if r.status_code == 200:
            time.sleep(PAUSE)
            return r.json()
        if r.status_code in (429, 500, 502, 503, 504):
            time.sleep(2 ** attempt * PAUSE)
            continue
        raise RuntimeError(f"STRING {endpoint} failed ({r.status_code}): {r.text[:300]}")
    raise RuntimeError(f"STRING {endpoint} failed after {retries} retries")


def _chunks(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def resolve(names, species=9606):
    """Map gene symbols / UniProt / Ensembl IDs to STRING.

    Returns {query: {"string_id": ..., "name": preferred_name}}; unmatched queries are absent.
    """
    out = {}
    for chunk in _chunks(list(names), 200):
        rows = _post("get_string_ids", {
            "identifiers": "\r".join(chunk), "species": species,
            "limit": 1, "echo_query": 1})
        for row in rows:
            out.setdefault(row["queryItem"],
                           {"string_id": row["stringId"], "name": row["preferredName"]})
    return out


def interaction_partners(string_ids, species=9606, required_score=700,
                         network_type="physical", limit=1000, min_escore=0.0):
    """All partners of each protein (partners do NOT need to be in the input list).

    required_score: combined score 0-1000 (700 = high confidence).
    min_escore: keep only edges whose *experimental* channel score is >= this (0-1).
    """
    edges = []
    for chunk in _chunks(list(string_ids), 50):
        rows = _post("interaction_partners", {
            "identifiers": "\r".join(chunk), "species": species,
            "required_score": required_score, "network_type": network_type,
            "limit": limit})
        for r in rows:
            if (r.get("escore") or 0) < min_escore:
                continue
            edges.append({
                "a": r["preferredName_A"], "b": r["preferredName_B"],
                "id_a": r["stringId_A"], "id_b": r["stringId_B"],
                "source": "STRING", "weight": float(r["score"]),
                **{k: r.get(k) for k in EVIDENCE}})
    return edges
