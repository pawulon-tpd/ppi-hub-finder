# ppi-hub-finder

Give it a list of proteins. It finds **which single protein, inside or outside your list,
directly interacts with the most of them**, and draws the interaction map.

Data comes from two public databases:

| Source | What it contains | Key needed |
|---|---|---|
| [STRING](https://string-db.org) | Combined evidence: experiments, curated complexes, text mining, co-expression, predictions. Scored 0–1. | No |
| [BioGRID](https://thebiogrid.org) | Experimentally curated interactions only (Y2H, AP-MS, co-IP…). | Free key |

## Install

```bash
git clone https://github.com/<you>/ppi-hub-finder.git
cd ppi-hub-finder
pip install -r requirements.txt
```

For BioGRID, register at <https://webservice.thebiogrid.org/> and:
```bash
export BIOGRID_ACCESS_KEY=your_key
```

## Usage

```bash
# STRING, physical network, high confidence (default)
python ppi_hubs.py TP53 MDM2 BAX BCL2 EGFR KRAS BRAF

# Names from a file, both databases, require experimental evidence in STRING
python ppi_hubs.py -f examples/genes_example.txt --source both --min-escore 0.4

# Only look at links between proteins in your own list
python ppi_hubs.py -f examples/genes_example.txt --only-within

# Also count proteins reached through one intermediate; best pair, all combinations tested
python ppi_hubs.py -f examples/genes_example.txt --depth 2 --cover 2 --exact

# Mouse
python ppi_hubs.py Trp53 Mdm2 Bax --species 10090
```

Inputs can be gene symbols, UniProt accessions or Ensembl IDs. STRING converts them
to canonical gene names, which BioGRID then uses too.

Main options:

| Option | Meaning |
|---|---|
| `--source string\|biogrid\|both` | Which database(s) to use |
| `--score 700` | STRING combined score cut-off (400 medium, 700 high, 900 highest) |
| `--network-type physical` | `physical` = binding/complex; `functional` = any association |
| `--min-escore 0.4` | Keep STRING edges only with experimental support |
| `--min-pubs 2` | BioGRID pair must be reported in ≥ N papers |
| `--only-within` | Ignore partners outside your list |
| `--depth 2` | Also count inputs reached through one intermediate protein |
| `--hop-weight 0.5` | Value of a 2-hop link relative to a direct link |
| `--max-shell 300` | Depth 2: how many first-shell partners are expanded |
| `--cover 3` | Combine 3 different targets that together reach the most of your list |
| `--exact` | Test every combination (among `--candidates` best singles) instead of greedy |

## Outputs (`results/`)

| File | Content |
|---|---|
| `target_ranking.csv` | Every protein ranked by direct links to your list: `n_input_partners`, `coverage_pct` (links + itself if in the list), and the names |
| `cover.csv` | Best combination of N different targets and what each one adds |
| `input_partner_counts.csv` | Number of interaction partners of each input protein, per database |
| `edges.csv` | All interactions with evidence scores |
| `network.png` | Static map: blue circles = your proteins, orange diamonds = outside proteins linking ≥2 of them, red = best target and its links, thick grey = found in both databases |
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
