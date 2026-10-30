# ppi-hub-finder

Give it a list of proteins. It finds **which protein, inside or outside your list,
interacts with the most of them** (directly or through one intermediate), picks the best
combination of N different targets, and draws the interaction map.

With `--mechanism` it also tells you **how** each target acts on your proteins (activates or
inhibits, by phosphorylation, ubiquitination, transcription, ...), directly or through one
intermediate protein.

Data comes from public databases:

| Source | What it contains | Key needed |
|---|---|---|
| [STRING](https://string-db.org) | Combined evidence: experiments, curated complexes, text mining, co-expression, predictions. Scored 0–1. | No |
| [BioGRID](https://thebiogrid.org) | Experimentally curated interactions only (Y2H, AP-MS, co-IP…). | Free key |
| [OmniPath](https://omnipathdb.org) | `--mechanism`: direction, activation/inhibition and modification sites, combined from ~100 resources (SIGNOR, PhosphoSitePlus, KEGG, Reactome, CollecTRI…). Human, mouse, rat. | No |
| [RCSB PDB](https://www.rcsb.org) + [UniProt](https://www.uniprot.org) | `--structures`: 3D structures, resolution and bound ligands per protein and per pair. `--ternary` adds ligand burial, exit vectors and clash-free poses. | No |

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

# ... and how each top target acts on your list (activates / inhibits, mechanism)
python ppi_hubs.py -f examples/genes_example.txt --source both --exclude-hubs --depth 2 --mechanism

# Druggability: PDB entries and ligands for the top targets and the best 3 pairs
python ppi_hubs.py -f examples/genes_example.txt --source both --exclude-hubs --structures

# ... plus ligand burial, exit vectors and clash-free ternary poses (PROTAC triage)
python ppi_hubs.py -f examples/genes_example.txt --source both --exclude-hubs --structures --ternary

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

**Mechanisms (OmniPath)**

| Option | Meaning | Default |
|---|---|---|
| `--mechanism` | Annotate how proteins act on each other; adds `paths.csv`, `target_effects.csv`, `mechanisms.csv` and colours the interactive map | off |
| `--mech-targets` | How many top targets get a path analysis and their own view in `network.html` | 20 |
| `--include-predicted` | Also keep records supported only by text mining or kinase-substrate prediction | off |
| `--omnipath-license` | `academic` or `commercial`: which OmniPath records you may use | academic |
| `--refresh-omnipath` | Re-download the cached OmniPath tables (done automatically after 30 days) | off |

**Structures (RCSB PDB)**

| Option | Meaning | Default |
|---|---|---|
| `--structures` | Look up PDB entries and bound ligands for the top targets and the best pairs | off |
| `--struct-targets` | How many top targets to look up | 10 |
| `--pairs` | How many best target pairs to report | 3 |
| `--ternary` | For each pair with a ligand on both sides: ligand burial, exit vectors, clash-free pose scan | off |
| `--ternary-vectors` | Exit vectors per ligand to try | 3 |
| `--ternary-poses` | Distinct arrangements written as PDB files per pair | 5 |
| `--refresh-structures` | Re-download cached PDB metadata (automatic after 30 days) | off |

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
| `network.html` | Interactive map (see below) |
| `paths.csv` | `--mechanism`: for each top target and each input it reaches, the route (direct or via a bridge), the effect on the input, the reverse effect, and each step's mechanism |
| `target_effects.csv` | `--mechanism`: per target, how many inputs it activates / inhibits / mixed / sign unknown / no data, with names |
| `mechanisms.csv` | `--mechanism`: every directed OmniPath record among the proteins in the network |

### Structural output (`--structures`, `--ternary`)

| File | Content |
|---|---|
| `top_pairs.csv` | The best `--pairs` target pairs by combined coverage (each protein appears in at most one pair) |
| `structures.csv` | Per protein: number of PDB entries, how many have a real ligand, best resolution, best ligand-bound entry, distinct ligands |
| `ligands.csv` | Every non-additive ligand found, with its PDB entry, chain, weight, resolution and name |
| `pair_structures.csv` | Per pair: the best ligand-bound entry on each side, and any entry containing **both** proteins |
| `ternary_summary.csv` | Per pair: ligand burial on each side, best attachment atoms and linker length, % of poses that are clash-free, number of distinct arrangements, warnings |
| `exit_vectors.csv` | Candidate linker attachment atoms: relative solvent exposure, SASA, how blocked the outward cone is, coordinates |
| `ternary_poses.csv` | The scanned poses (top 200 per pair) with linker length, tilt, spin, clashes and contacts |
| `poses/pose_<A>_<B>_<n>.pdb` | The distinct arrangements as PDB files — open in PyMOL or ChimeraX |

### The interactive map (`network.html`)

Open it in any browser (needs internet for the graph library).

- **View:** *Overview* shows your proteins and the best outside proteins. Pick a target from the
  list, or double-click a protein, to see **only its routes to your list**; the side panel then
  lists every input with its route(s) and the overall effect. Hover a row for each step's mechanism.
- **Effect toggles:** activates (green arrow), inhibits (red T-bar), mixed (purple),
  sign unknown (dark grey arrow), binding / association only (light grey line).
- **Mechanism toggles:** phosphorylation, ubiquitination, acetylation, transcriptional regulation,
  other regulation, ... (only those present in your network).
- **Links only in OmniPath:** regulatory links that STRING/BioGRID do not have (dashed, off by default).
- Hover any protein or link for details (sites, databases, number of references).

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

### Mechanisms (`--mechanism`)

For each top target *T* and each input *Y* it reaches:

- **Direct:** OmniPath records for *T → Y* (effect on your protein) and *Y → T* (reverse).
- **Via a bridge *X*:** the chain *T → X → Y*. Its effect is the product of both steps:
  activate × inhibit = inhibit, inhibit × inhibit = activate. Up to 3 bridges are kept per input,
  best first: complete signed chain, then clear sign over mixed, then the chain whose weakest
  step has the most references.
- A **direct partner with no directed record** also gets bridges, because binding does not say
  how one acts on the other (e.g. AKT1 binds TP53, but acts on it via AKT1 → MDM2 ⊣ TP53).
- If the routes to one input disagree, the overall effect is **mixed**. Look at `paths.csv`.

Records supported only by text mining or prediction (REACH, Sparser, RLIMS-P, MIMP, NetworKIN,
BEL Large Corpus) are dropped unless `--include-predicted`. Protein complexes are skipped.
A chained effect is a **hypothesis**: it depends on cell type, timing and which sites are hit.

The two OmniPath tables (~40 MB for human) are downloaded once and cached in
`~/.cache/ppi-hub-finder/`.

### How the structural part works

**Which ligands count.** Crystallisation additives, buffers, cryoprotectants, glycans and
ions (`structures.ADDITIVES`) are never counted as ligands. Cofactors and nucleotides
(ATP, NAD, HEM…) are tracked separately, because an ATP site is a real pocket but not a
drug-like starting point on its own. Anything else is `fragment` (150–300 Da) or
`drug-like` (> 300 Da). **A ligand only counts for the protein whose chains it sits on**,
so in a co-complex the partner is not credited with it.

**Ligand burial** is Shrake-Rupley SASA (Biopython): the ligand's exposed area in the
complex divided by its area free in solution. 0.99 means a fully enclosed pocket with
nowhere for a linker to leave; ~0.6–0.9 is typical of a pocket with a solvent-facing rim.

**Exit vectors** are the ligand's solvent-exposed heavy atoms, each with a direction
pointing away from the local protein wall, and an `openness` score counting protein atoms
inside a 35° cone, 12 Å long, along that direction. If no atom is exposed, the most
exposed atoms are still returned but flagged `no_open_vector` — that is itself the answer.

**Ternary pose scan.** Protein A is fixed. Protein B is placed so its exit vector points
back at A's, with four parameters sampled: the attachment-point separation (the linker
span, 5–20 Å), the tilt of the A→B direction away from A's exit vector (0/20/40°), the
azimuth of that tilt, and B's spin about the axis. Each pose is scored by heavy-atom
clashes (< 2.5 Å) and interface contacts (< 4.5 Å). Thousands of near-identical poses are
then clustered so you get a handful of genuinely different arrangements.

### How much to trust the ternary scan

`validation/benchmark_5t35.py` runs the scan on the two binary structures (BRD4+JQ1,
VHL+VH032) and compares every pose with the real BRD4:MZ1:VHL ternary complex (5T35):

```
poses 16848 | clash-free 11917 (71%) | distinct arrangements 10
sampling  - closest pose to native:        5.4 A
scoring   - top-ranked pose:               30.6 A
            median over clash-free poses:  30.8 A
```

Read that honestly:

- **The sampling works.** The native arrangement is inside the searched set (5.4 Å CA RMSD).
- **The scoring does not.** Ranking by clashes and contacts does not pick it out; the
  top-ranked pose is 30 Å away. Clash-freedom is weak evidence — 71% of poses pass.
- **So this is triage, not prediction.** A pair where *no* linker length gives a clash-free
  pose, or where a ligand is ~99% buried, is a real red flag. A pair with many clash-free
  arrangements is permissive — it does not tell you the PROTAC will work, which depends on
  linker chemistry, induced fit and cooperativity that no rigid scan captures.
- Everything is computed from **one crystal conformation** of each protein, and proteins
  move. Use the pose files as hypotheses to model properly, not as answers.

## Caveats

- STRING's score is not pure experiment. Use `--min-escore` or BioGRID for "proven binding".
- High-throughput screens make well-studied proteins look like hubs (study bias). Check
  `n_publications` / `systems` in `edges.csv`.
- Outside proteins only show links to your list, not their full partner count.
- Partners per protein are capped by `--limit` (STRING).
- An interaction is not druggability. The best hub may be essential or undruggable.

## Data licences

RCSB PDB data is in the public domain (CC0); UniProt is CC BY 4.0.
STRING: CC BY 4.0. BioGRID: MIT. OmniPath: each record keeps the licence of its original
resource; `--omnipath-license commercial` keeps only records allowed for commercial use.
Please cite them if you publish results.

## Tests

```bash
pip install pytest && pytest
```
