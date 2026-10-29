# Handoff: ppi-hub-finder

## Goal
Given a list of proteins (gene symbols / UniProt / Ensembl), find which protein, in the list
or outside it, interacts with the most of them (directly or via one intermediate), pick the
best combination of N different targets, and draw the interaction map. The owner (Paweł)
does drug-discovery work. This is a GitHub portfolio repo, so keep it clean, documented
and tested.

## Status
- Written in a claude.ai chat sandbox without network access, then moved here.
- 2026-10-29: first live STRING run works (depth 1 and depth 2 + `--exact`, ~15 s for
  the 14-gene example). `get_string_ids` echoes `queryItem`, `score` is 0-1 as assumed.
- **BioGRID live calls still untested** (needs `BIOGRID_ACCESS_KEY`).
- Tests: `pytest` (config in `pyproject.toml` puts the repo root on `sys.path`).

## Layout
- `ppi_hubs.py`: CLI entry point (argparse). Flow: resolve names (STRING) → fetch
  partners → optional depth-2 shell expansion → dedupe → graph → rank → cover → plots → CSVs.
- `ppi/string_api.py`: STRING REST (`get_string_ids`, `interaction_partners`), POST,
  1 s pause, `caller_identity`, retries. Edges keep `id_a/id_b` (needed for depth 2).
- `ppi/biogrid_api.py`: BioGRID REST, key from `BIOGRID_ACCESS_KEY`, paging (max 10000),
  dedupe by interaction id, physical-only filter, collapses records to one edge per pair
  (weight = number of distinct PubMed IDs).
- `ppi/network.py`: `dedupe_edges`, `build_graph` (edge attr `sources` = set of DBs),
  `partner_counts`, `reach_table` ({protein: {input: hops}}, hops 0 self / 1 direct / 2 via
  bridge), `rank_targets` (score = direct + hop_weight × two_hop), `greedy_cover`,
  `exact_cover` (all combos among top-N candidates; refuses > 5M combos), `two_hop_paths`.
- `ppi/plot.py`: matplotlib PNG and pyvis HTML. Blue circles = inputs, orange diamonds =
  outside proteins, red = best target (solid = direct, dashed = 2-hop paths), thick grey =
  edge in both DBs.
- `tests/test_network.py`, `examples/genes_example.txt`, `README.md` (user docs).

## Key design decisions
- Cover targets may be any protein in the graph, are always distinct, and an input reached
  by several targets counts once at its best distance (weighted max coverage).
- Depth 2 only expands the top `--max-shell` (300) first-shell partners, so targets are
  limited to inputs + first shell. Second-shell targets are deliberately not considered (cost).
- STRING names are canonical; BioGRID names are mapped case-insensitively onto them.

## Things to verify on first live run
1. `python ppi_hubs.py -f examples/genes_example.txt`: does `get_string_ids` with
   `echo_query=1` return `queryItem`? Is `score` 0-1 in JSON? Are `escore`, etc. present?
2. `--source both`: BioGRID response shape (dict keyed by id), field names
   `OFFICIAL_SYMBOL_A/B`, `EXPERIMENTAL_SYSTEM_TYPE`, `PUBMED_ID`; does paging work?
3. `--depth 2`: runtime and memory with 300 shell proteins × limit 1000; check the
   pyvis tooltips (`\n` may render literally, in which case switch to `<br>`).

## Ideas / possible next steps
- Option to exclude promiscuous hubs (UBC, HSP90AA1, ...) or penalise by total degree.
- Cache API responses (e.g. `requests-cache`) to avoid re-fetching.
- Edge-confidence-weighted paths at depth 2 (product of scores).
- Add structural evidence (PDBe/RCSB complexes) as a third source.
- GitHub Actions running pytest.

## Commands
pip install -r requirements.txt pytest
pytest -q
python ppi_hubs.py -f examples/genes_example.txt --depth 2 --cover 2 --exact
