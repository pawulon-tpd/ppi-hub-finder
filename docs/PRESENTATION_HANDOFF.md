# Handoff: building the ppi-hub-finder slide deck

Everything needed to assemble the PowerPoint on another machine. Figures are in this
folder as PNGs at 200 dpi, sized for 16:9. Regenerate any of them with:

```bash
python docs/make_figures.py
```

| Figure | Use it for |
|---|---|
| `pipeline.png` | The one-slide overview of the whole tool |
| `ternary_scan.png` | Explaining how the PROTAC geometry search works |
| `benchmark_5t35.png` | The validation / honesty slide |
| `results/*/network.png` | A real interaction map (generate from your own gene list) |

---

## 1. The one-sentence version

> Give it a list of proteins. It finds which protein — inside or outside the list —
> interacts with the most of them, works out how it acts on them, checks what structures
> and ligands exist, and asks whether a degrader could physically link the two.

The thread running through the talk: **each stage narrows a list.** Hundreds of
interaction partners become a handful of ranked targets, then a few with a real chemical
starting point, then the ones where the geometry actually permits a ternary complex.

---

## 2. The pipeline, stage by stage

Use `pipeline.png` while talking through this. Stages 1–3 always run; 4, 5 and 6 are
opt-in flags.

### Stage 1 — Resolve names (STRING)
Input can be gene symbols, UniProt accessions or Ensembl IDs. STRING maps them to
canonical gene names, which everything downstream then shares.

*Why it matters:* BioGRID still uses some outdated symbols — it calls parkin **PARK2**
where STRING says **PRKN**. Without reconciling these, one protein becomes two nodes and
its links are split in half. The tool matches them through BioGRID's own synonym lists.

### Stage 2 — Fetch interactions (STRING + BioGRID)
- **STRING**: combined evidence (experiments, curated complexes, text mining,
  co-expression). Scored 0–1, no key needed.
- **BioGRID**: experimentally curated interactions only. Free key.
- Optional **depth 2**: also count proteins reached through one intermediate.

*Two decisions worth defending on a slide:*
- **`--min-pubs 2` is the default.** Requiring two independent publications stops
  single high-throughput screens from dominating. With it off, PRKN, TRIM67 and NUDT21
  top the list; with it on, the biology-driven ranking returns.
- **`--exclude-hubs`** removes 23 proteins that stick to everything — HSP90/HSP70
  chaperones, ubiquitin, 14-3-3. They are real biology but they drown out the specific
  links you are looking for.

### Stage 3 — Rank targets and pairs
Each protein is scored by how much of *your* list it reaches: direct partners, plus
2-hop partners at half weight. Pairs are then chosen so that together they cover the
most of your list, each protein used at most once.

*Key idea for the slide:* coverage is **weighted max-coverage** — an input reached by
both members of a pair counts once, at its best distance. So a pair is only as good as
what the second member *adds*.

### Stage 4 — Mechanism (OmniPath) · `--mechanism`
Interaction databases say *that* two proteins touch. OmniPath says *what one does to the
other*: direction, activation or inhibition, and the modification with its exact residue.

For two-step routes the effects multiply — inhibit × inhibit = activate. The textbook
example comes straight out of the data:

```
AKT1 --phosphorylates S166, activates--> MDM2 --ubiquitinates K373, inhibits--> TP53
                        so AKT1 indirectly inhibits TP53
```

*Caveats to state out loud:* many links have no direction at all; some are flagged both
activating and inhibiting (context-dependent); and records supported only by automated
text mining are dropped by default, because they contain errors like "MDM2
phosphorylates TP53" — MDM2 is not a kinase.

### Stage 5 — Structures and ligands (RCSB PDB + UniProt) · `--structures`
Per protein: how many PDB entries exist, how many have a real ligand, the best
resolution, the distinct ligands. Per pair: the best ligand-bound structure on each side,
and whether any structure contains **both** proteins.

*Two things the tool is careful about, both worth a sentence:*
- Crystallisation additives, buffers and ions are never counted as ligands; cofactors
  (ATP, NAD, HEM) are tracked separately, because an ATP site is a real pocket but not a
  drug-like starting point on its own.
- **A ligand only counts for the protein whose chains it sits on.** In a co-complex the
  partner would otherwise be credited with a ligand it does not bind. Fixing this changed
  real output: TP53's ligand-bearing structure count dropped from 95 to 65, and BRAF's
  best structure moved from one whose ligand is 99.6% buried and unusable to one at 62%.

### Stage 6 — Ternary geometry · `--ternary`, `--e3`
Use `ternary_scan.png` here. For a pair with a ligand on each side:

1. **Ligand burial** — Shrake-Rupley solvent-accessible surface, the ligand in the
   complex versus free in solution. 99% buried means a fully enclosed pocket with
   nowhere for a linker to leave.
2. **Exit vectors** — the solvent-exposed ligand atoms, each with a direction pointing
   away from the protein wall and a score for how blocked that direction is.
3. **Pose scan** — protein A fixed, protein B placed so the exit vectors face each
   other; sweep linker span (5–20 Å), tilt, azimuth and spin. Score every pose by
   heavy-atom clashes and contacts, then cluster away near-duplicates.

Output: what fraction of poses are clash-free, and a handful of genuinely different
arrangements written as PDB files you can open in PyMOL.

### E3 ligase mode · `--e3`
Instead of only pairing proteins from your list, pair each top target with E3 ligases
that have published small-molecule handles: **CRBN, VHL** (clinical), **MDM2, XIAP,
BIRC2** (established), **KEAP1** and others (emerging). Only gene symbols are hard-coded —
structures and ligands are looked up live, so the table does not go stale.

Each E3 carries its on-target caveat, which is easy to forget mid-analysis: IMiD binding
to CRBN is poor in mouse and recruits neosubstrates; IAP antagonists auto-degrade; MDM2
binders stabilise p53; KEAP1 binders stabilise NRF2.

---

## 3. Numbers you can quote

From the 14-gene example (apoptosis + MAPK/PI3K), STRING + BioGRID, hubs excluded:

| | |
|---|---|
| Network size | 1 841 proteins, 4 988 interactions from 14 inputs |
| Best single target | TP53 — 10 of 14 direct, 3 more via one intermediate |
| Best pair | TP53 + BRAF — 100% of the list |
| Runtime, 24 proteins | 25 s (direct) / 90 s (with 2-hop) |
| First structural run | ~2 000 PDB entries, ~6 min, cached after |
| Ternary scan | ~17 000 poses per vector pair, 36 s → 9 s on 10 cores |

PROTAC tractability for TP53, from `e3_pairs.csv`:

| E3 | Target ligand | E3 ligand | Clash-free poses |
|---|---|---|---|
| VHL | 6GGC / EXN, 92% buried | 9GIO / 3JF, 86% | 81.7% |
| CRBN | 6GGC / EXN | 9CUO, 75% | 79.8% |
| MDM2 | 6GGC / EXN | 6Q9L / HTZ, 86% | 68.2% |

A nice sanity check to mention: TP53–MDM2 came back with 2 co-complex structures, which
is correct — those are the real p53:MDM2 crystals.

---

## 4. The honesty slide — do not skip it

Use `benchmark_5t35.png`. This is the slide that earns trust in the rest.

The only fair way to judge a pose scan is against a crystal structure of the thing it
approximates. **5T35** is BRD4 + the PROTAC MZ1 + VHL. The test takes the two *binary*
structures the scan would normally start from (BRD4+JQ1, VHL+VH032), runs the scan, and
measures every pose against the real complex.

```
poses 16 848 | clash-free 11 917 (71%)
sampling - closest pose to native:   5.4 A
scoring  - top-ranked pose:         30.6 A
           median clash-free pose:  30.8 A
```

Two separate conclusions, and the distinction is the whole point:

- **The sampling works.** The native arrangement is inside the searched set.
- **The scoring does not.** Clashes and contacts do not single it out — and with 71% of
  poses clash-free, clash-freedom is weak evidence on its own.

So the honest claim is: **a pair where no linker length gives a clash-free pose, or where
a ligand is ~99% buried, is a real red flag. A pair with many clash-free arrangements is
permissive — not a prediction that the PROTAC will work.** Everything is computed from a
single crystal conformation, and proteins move.

Reproduce it live with `python validation/benchmark_5t35.py`.

---

## 5. Suggested slide order (~14 slides)

| # | Slide | Visual |
|---|---|---|
| 1 | Title | — |
| 2 | The problem: a list of proteins, which one do I go after? | — |
| 3 | What the tool does, in one sentence | — |
| 4 | Pipeline overview | `pipeline.png` |
| 5 | Data sources and what each is good for | table from §2 |
| 6 | Ranking: direct + 2-hop coverage, best pairs | `results/*/network.png` |
| 7 | Two decisions that change the answer: `--min-pubs 2`, `--exclude-hubs` | before/after table |
| 8 | Mechanism: not *that* they interact but *what* one does | the AKT1→MDM2⊣TP53 chain |
| 9 | Structures and ligands: is there a chemical starting point? | `structures.csv` extract |
| 10 | Ternary geometry: how the scan works | `ternary_scan.png` |
| 11 | E3 mode: which target, which ligase | `e3_pairs.csv` extract |
| 12 | **Validation and limits** | `benchmark_5t35.png` |
| 13 | What it is good for / what it is not | — |
| 14 | Next step: a real scoring function | — |

---

## 6. Likely questions, and honest answers

**"Why not just use STRING's website?"** STRING ranks interactions; this ranks *targets
against your list*, combines several targets into a covering set, adds mechanism and
direction, and carries it through to structures and geometry.

**"Is the 2-hop step real biology?"** Treat it as "possible indirect influence", not
binding. Interactomes are small worlds — at two steps, big hubs reach almost everything,
which is exactly why `--exclude-hubs` and a high confidence cut-off matter.

**"Could it find a protein that doesn't touch my list directly?"** No, deliberately. Both
the target and the intermediate must be in the first shell. Going further would mean
downloading the partners of ~1 800 proteins, and at that distance everything connects to
everything.

**"Can it predict whether my PROTAC will work?"** No. It can tell you whether a linker
could physically span the two ligands without the proteins overlapping. Linker chemistry,
induced fit and cooperativity are not modelled — see the validation slide.

**"Can I use this commercially?"** STRING (CC BY 4.0), BioGRID (MIT), RCSB PDB (public
domain) and UniProt (CC BY 4.0) are all fine with citation. **OmniPath needs care** —
it aggregates ~230 resources, each keeping its own licence, and several large ones
(PhosphoSitePlus, HPRD, phospho.ELM) are non-commercial. The tool defaults to
`--omnipath-license commercial`, which keeps only records a commercially-licensed
resource supports. Note this filter runs in the tool, not on the server: OmniPath's own
`license` parameter was tested and returned identical data either way.

---

## 7. Commands to produce fresh figures for the deck

```bash
cd ppi-hub-finder && source .venv/bin/activate
set -a; . ./.env; set +a          # BioGRID key
```

A network map plus the full structural output from your own list:

```bash
python ppi_hubs.py -f my_targets.txt --source both --exclude-hubs --depth 2 \
  --mechanism --structures --e3 --ternary --out results/talk
```

That writes `results/talk/network.png` (static map), `network.html` (interactive, good
for a live demo), the CSV tables quoted above, and `poses/*.pdb` for a PyMOL screenshot.

For a PyMOL figure of a predicted arrangement:

```
load results/talk/poses/pose_e3_TP53_VHL_1.pdb
hide everything; show cartoon; show sticks, hetatm and not resn HOH
color grey80; color skyblue, chain A; color orange, chain b
bg_color white; set ray_opaque_background, 1; orient; ray 2000, 1400
png ternary_pose.png, dpi=200
```
