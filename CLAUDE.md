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
- BioGRID live run works too (key in gitignored `.env`; `set -a; . ./.env; set +a`).
  Field names correct, paging needed and works (~24k records for the example).
  BioGRID keeps some outdated symbols (PARK2 vs STRING's PRKN): `biogrid_api._namer`
  maps them onto STRING names via `SYNONYMS_A/B` (only if exactly one synonym matches).
- `--min-pubs` defaults to 2: with 1, BioGRID high-throughput screens dominate the ranking
  (PRKN, TRIM67, NUDT21 on top).
- `--exclude-hubs [NAMES]` drops `network.PROMISCUOUS` (chaperones, ubiquitin, 14-3-3)
  plus extra names from the edges (never inputs), so they are neither targets nor bridges.
- Depth 2 + BioGRID asks only for links among inputs + shell (`includeInteractors=false`,
  one request): ~50 s for the example. Fetching all shell partners took >10 min.
- Tests: `pytest` (config in `pyproject.toml` puts the repo root on `sys.path`).

- 2026-10-30: `--mechanism` (OmniPath) + new interactive `network.html` (custom vis-network
  page with toggles, replaced pyvis). Tested live on the example (human, depth 1/2), mouse,
  `--only-within`, unsupported species; page checked in a browser (no console errors).

- 2026-10-30 (later): `--structures` (RCSB/UniProt inventory) and `--ternary` (ligand SASA,
  exit vectors, clash-checked rigid pose scan) + `validation/benchmark_5t35.py`.
  **Local only, not pushed** - Pawel wants to try it on a real target list first.

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
- `ppi/mechanisms.py`: OmniPath `interactions` (datasets omnipath + collectri) and `enz_sub`
  tables, cached in `~/.cache/ppi-hub-finder` (30 days). `build_index` -> {(src, tgt): sign,
  kinds, details}. Sign: consensus_* first, else is_*; several records disagree -> mixed (0).
  Drops text-mined/predicted-only records (`LOW_CONFIDENCE`) and COMPLEX rows.
  `target_paths`: per target/input, direct + bridge chains (sign = product), bridges ranked
  by `_bridge_rank`; direct partners without a record also get bridge chains.
- `ppi/interactive.py`: builds JSON + self-contained vis-network page (CDN). Overview vs
  per-target view (DataView filters), effect/mechanism toggles, side panel per target.
- `ppi/structures.py`: UniProt gene->accession, RCSB search (accession->entries) and
  GraphQL (resolution, chains, ligands), all cached per entry in `~/.cache/ppi-hub-finder/rcsb`.
  `ADDITIVES`/`COFACTORS` classify ligands; `ligands_of(entry, acc)` attributes a ligand to a
  protein by auth-chain intersection (without this, co-complexes credit both partners).
- `ppi/geometry.py`: mmCIF fetch+cache, `ligand_sasa` (Biopython Shrake-Rupley, complex vs
  free ligand), `exit_vectors` (exposed atoms + outward direction + cone openness; falls back
  to the most exposed atoms flagged `buried`), `ternary_scan` (samples linker length, tilt,
  azimuth, spin; counts clashes/contacts), `cluster_poses`, `write_pose`.
  `chains_near_ligand` drops crystal-packing neighbours by a 15% contact threshold.
  `_apply` wraps matmul in errstate: numpy 2.0 + Apple Accelerate emits spurious warnings.
- `ppi/plot.py`: matplotlib PNG. Blue circles = inputs, orange diamonds =
  outside proteins, red = best target (solid = direct, dashed = 2-hop paths), thick grey =
  edge in both DBs.
- `tests/test_network.py`, `examples/genes_example.txt`, `README.md` (user docs).

## Key design decisions
- Cover targets may be any protein in the graph, are always distinct, and an input reached
  by several targets counts once at its best distance (weighted max coverage).
- Depth 2 only expands the top `--max-shell` (300) first-shell partners, so targets are
  limited to inputs + first shell. Second-shell targets are deliberately not considered (cost).
- STRING names are canonical; BioGRID names are mapped case-insensitively onto them,
  falling back to BioGRID synonyms for outdated symbols.

## Structural work: what the benchmark says
`python validation/benchmark_5t35.py` (BRD4+JQ1 / VHL+VH032 vs the real 5T35 ternary):
sampling reaches the native pose (5.4 A) but scoring by clashes+contacts does NOT rank it
(top pose 30.6 A; 71% of poses are clash-free). So the output is triage, not prediction -
keep saying so in docs and output. Improving this means a real scoring function
(shape complementarity / desolvation / statistical potential), not more sampling.

## OmniPath licences
The server's `license=academic|commercial` parameter is ACCEPTED BUT IGNORED (verified
2026-10-30: identical row counts and PhosphoSite present in both). Filtering is therefore
done client-side in `mechanisms.commercial_ok` using `/resources` metadata (cached as
omnipath_resources.json): a record survives only if some source resource has
purpose == "commercial". Composite aggregators (CollecTRI, DoRothEA) and unlisted
resources do not count. Keeps ~92% of mechanisms on the 24-gene example.

## Open questions
- Real-world target lists not tested yet (only the 14-gene example).
- No E3-ligase mode yet: `--structures` ranks whatever is in the list. A built-in E3 table
  (VHL, CRBN, MDM2, XIAP, DCAF15...) with known ligand-bound PDBs would answer "which of my
  proteins is PROTAC-tractable, and with which E3" - Pawel raised this, not yet decided.
- First `--structures` run on ~20 proteins fetches ~2000 PDB entries (~6 min); cached after.
- BioGRID used to time out above ~10 proteins: partner pulls now go out in chunks of 10,
  within-shell pulls in chunks of 150 PLUS every chunk pair (else cross-chunk pairs are
  lost), PAGE 5000, and `_get` retries with a longer timeout each round. 24 proteins:
  25 s at depth 1, 90 s at depth 2 (`examples/genes_24.txt`).
- Overview in network.html gets dense with many inputs; per-target views are the readable part.

## Ideas / possible next steps
- Penalise hubs by total degree instead of (or besides) the fixed `--exclude-hubs` list.
- Cache API responses (e.g. `requests-cache`) to avoid re-fetching.
- Edge-confidence-weighted paths at depth 2 (product of scores).
- Option to use OmniPath directed links as a ranking source (now display/annotation only).
- Map OmniPath complexes (COMPLEX:...) onto member proteins.
- 2-hop paths of outside targets only reach inputs via first-shell bridges; second shell skipped.
- Add structural evidence (PDBe/RCSB complexes) as a third source.
- GitHub Actions running pytest.

## Commands
pip install -r requirements.txt pytest
pytest -q
python ppi_hubs.py -f examples/genes_example.txt --depth 2 --cover 2 --exact
set -a; . ./.env; set +a   # BioGRID key
python ppi_hubs.py -f examples/genes_example.txt --source both --exclude-hubs --depth 2 --mechanism
