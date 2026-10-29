# ppi-hub-finder

Give it a list of proteins. It finds **which protein, inside or outside your list,
interacts with the most of them** (directly or through one intermediate), picks the best
combination of N different targets, and draws the interaction map.

Data comes from two public databases:

| Source | What it contains | Key needed |
|---|---|---|
| [STRING](https://string-db.org) | Combined evidence: experiments, curated complexes, text mining, co-expression, predictions. Scored 0–1. | No |
| [BioGRID](https://thebiogrid.org) | Experimentally curated interactions only (Y2H, AP-MS, co-IP…). | Free key |

## Install

Needs Python 3.9+.

```bash
git clone https://github.com/pawulon-tpd/ppi-hub-finder.git
cd ppi-hub-finder
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

STRING needs no key. For BioGRID, register for a free key at
<https://webservice.thebiogrid.org/> and save it in `.env` (git-ignored):

```bash
echo 'BIOGRID_ACCESS_KEY=your_key_here' > .env
```

## Quick start

In every new terminal:

```bash
cd ppi-hub-finder && source .venv/bin/activate
set -a; . ./.env; set +a      # only needed for BioGRID
```

```bash
# A few genes typed directly (STRING only)
python ppi_hubs.py TP53 MDM2 BAX BCL2 EGFR KRAS

# Genes from a file (spaces, commas or new lines; # starts a comment)
python ppi_hubs.py -f examples/genes_example.txt

# Recommended everyday run: both databases, sticky hubs removed
python ppi_hubs.py -f examples/genes_example.txt --source both --exclude-hubs

# Full analysis: also count proteins reached via one intermediate, test every pair of targets
python ppi_hubs.py -f examples/genes_example.txt --source both --exclude-hubs --depth 2 --cover 2 --exact

# Mouse
python ppi_hubs.py Trp53 Mdm2 Bax --species 10090
```

Inputs can be gene symbols, UniProt accessions or Ensembl IDs. STRING converts them
to current gene names, which BioGRID then uses too.

Tip: use `--out results/some_name` for each run so earlier results are not overwritten.

## Options (`python ppi_hubs.py -h`)

**Input**

| Option | Meaning | Default |
|---|---|---|
| `proteins` | Gene names, UniProt or Ensembl IDs | – |
| `-f FILE` | Read names from a file | – |
| `--species` | NCBI taxon id: 9606 human, 10090 mouse, 10116 rat, 4932 yeast | 9606 |
| `--biogrid-taxid` | BioGRID taxon id if it differs (yeast: 559292) | same as `--species` |

**Databases and evidence**

| Option | Meaning | Default |
|---|---|---|
| `--source` | `string`, `biogrid` or `both` | string |
| `--score` | STRING confidence 0–1000 (400 medium, 700 high, 900 highest) | 700 |
| `--network-type` | `physical` = binding/complex; `functional` = any association | physical |
| `--min-escore` | Keep STRING links with experimental score ≥ this (0–1), e.g. 0.4 | 0 (off) |
| `--limit` | Max STRING partners per protein | 1000 |
| `--min-pubs` | BioGRID pair must be reported in ≥ N papers (1 also keeps single high-throughput screens, which tend to flood the ranking) | 2 |
| `--all-biogrid` | Also keep BioGRID genetic (non-binding) interactions | off |
| `--biogrid-key` | BioGRID key, if not set in the environment | `BIOGRID_ACCESS_KEY` |

**Filtering**

| Option | Meaning |
|---|---|
| `--exclude-hubs` | Drop sticky proteins from the network (list below). Your own input proteins are never dropped |
| `--exclude-hubs VCP SQSTM1` | The same, plus the names you add |
| `--only-within` | Only links between your own proteins; ignore all others |

Proteins removed by `--exclude-hubs` (human gene names, matched case-insensitively):

- **Chaperones:** HSP90AA1, HSP90AB1, HSP90B1, HSPA1A, HSPA1B, HSPA4, HSPA5, HSPA8, HSPA9, HSPD1, CDC37, STUB1
- **Ubiquitin:** UBB, UBC, UBA52, RPS27A
- **14-3-3:** YWHAB, YWHAE, YWHAG, YWHAH, YWHAQ, YWHAZ, SFN

Some of their links are real biology (HSP90 stabilises EGFR, AKT1 and BRAF; 14-3-3 regulates RAF1),
so compare runs with and without the flag.

**Search**

| Option | Meaning | Default |
|---|---|---|
| `--depth 2` | Also count inputs reached through one intermediate protein | 1 |
| `--hop-weight` | Value of a 2-hop link relative to a direct one | 0.5 |
| `--max-shell` | Depth 2: how many direct partners are explored further | 300 |
| `--cover N` | Find the best combination of N different targets | 3 |
| `--exact` | Test every combination (among `--candidates` best singles) instead of greedy | off |
| `--candidates` | `--exact`: how many top single targets to combine | 60 |

**Output**

| Option | Meaning | Default |
|---|---|---|
| `--out DIR` | Output folder | `results` |
| `--top N` | Rows printed on screen | 10 |
| `--max-external` | Max outside proteins drawn on the map | 25 |

## Outputs

| File | Content |
|---|---|
| `target_ranking.csv` | Every protein ranked by how much of your list it reaches: `n_input_partners`, `n_two_hop`, `score`, `coverage_pct` and the names |
| `cover.csv` | Best combination of N different targets and what each one adds |
| `input_partner_counts.csv` | Number of interaction partners of each input protein, per database |
| `edges.csv` | All interactions with evidence scores |
| `network.png` | Static map: blue circles = your proteins, orange diamonds = outside proteins linking ≥2 of them, red = best target and its links (dashed = 2-hop paths), thick grey = found in both databases |
| `network.html` | Interactive version (drag, zoom, hover for details) |

## How the ranking works

For every protein *P* in the fetched network:

- **n_input_partners**: proteins from your list that *P* directly interacts with
- **n_two_hop** (depth 2): inputs reached only through one intermediate (*P → X → input*);
  `two_hop_partners` names the intermediates
- **score** = n_input_partners + hop_weight × n_two_hop
- **coverage**: everything reached, plus *P* itself if it is in your list

### Combining targets (`--cover N`)

Targets can be any protein in the network, from your list or outside it. They are always
different proteins. An input reached by several targets counts once, at its best distance,
so overlap is never double-counted. `newly_reached` shows what each target adds.

- **Greedy (default):** best single target, then the one adding the most not-yet-reached
  inputs, and so on. Fast and usually near-optimal, but it can miss a better pair.
- **`--exact`:** tests every combination of N among the top `--candidates` single targets.

### Depth 2

After fetching the inputs' partners (first shell), the tool fetches interactions among the
top `--max-shell` first-shell proteins. This lets an outside protein reach an input through
a bridge. Proteins outside the first shell are not considered as targets.

Interactomes are "small worlds": at 2 hops, big hubs (ubiquitin, chaperones) reach almost
everything. Keep a high `--score`, use `--min-escore` or BioGRID, and read 2-hop hits as
"possible indirect influence", not binding.

## Caveats

- STRING's score is not pure experiment. Use `--min-escore` or BioGRID for "proven binding".
- High-throughput screens make well-studied proteins look like hubs (study bias). Check
  `n_publications` / `systems` in `edges.csv`.
- Outside proteins only show links to your list, not their full partner count.
- Partners per protein are capped by `--limit` (STRING).
- An interaction is not druggability. The best hub may be essential or undruggable.

## Data licences

STRING: CC BY 4.0. BioGRID: MIT. Please cite them if you publish results.

## Tests

```bash
pip install pytest && pytest
```
