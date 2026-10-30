#!/usr/bin/env python3
"""Give a list of proteins -> find which protein (in or outside the list)
reaches the most of them, directly or through one intermediate, and draw the map.

Examples:
  python ppi_hubs.py TP53 MDM2 BAX BCL2 EGFR KRAS
  python ppi_hubs.py -f examples/genes_example.txt --source both --min-escore 0.4
  python ppi_hubs.py -f examples/genes_example.txt --depth 2 --cover 2 --exact
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

from ppi import biogrid_api, interactive, mechanisms, network, plot, string_api, structures


def parse():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("proteins", nargs="*", help="gene names / UniProt / Ensembl IDs")
    p.add_argument("-f", "--file", help="text file with names (spaces, commas or new lines; # comments)")
    p.add_argument("--source", choices=["string", "biogrid", "both"], default="string")
    p.add_argument("--species", type=int, default=9606,
                   help="NCBI taxon id: 9606 human, 10090 mouse, 10116 rat, 4932 yeast (STRING)")
    p.add_argument("--biogrid-taxid", type=int, help="BioGRID taxon id if different (yeast: 559292)")
    p.add_argument("--score", type=int, default=700, help="STRING min combined score 0-1000 [700]")
    p.add_argument("--network-type", choices=["physical", "functional"], default="physical")
    p.add_argument("--min-escore", type=float, default=0.0,
                   help="STRING min experimental-channel score 0-1, e.g. 0.4 [0 = off]")
    p.add_argument("--limit", type=int, default=1000, help="STRING max partners per protein [1000]")
    p.add_argument("--biogrid-key", help="BioGRID access key (or env BIOGRID_ACCESS_KEY)")
    p.add_argument("--min-pubs", type=int, default=2,
                   help="BioGRID min publications per pair; 1 also keeps single high-throughput "
                        "screens [2]")
    p.add_argument("--all-biogrid", action="store_true", help="also keep BioGRID genetic interactions")
    p.add_argument("--exclude-hubs", nargs="*", metavar="NAME",
                   help="drop sticky proteins (chaperones, ubiquitin, 14-3-3) from the network; "
                        "extra names are added to that list. Your input proteins are never dropped")
    p.add_argument("--only-within", action="store_true",
                   help="ignore partners outside your list (only rank your own proteins)")
    p.add_argument("--depth", type=int, choices=[1, 2], default=1,
                   help="1 = direct partners; 2 = also inputs reached via one intermediate [1]")
    p.add_argument("--hop-weight", type=float, default=0.5,
                   help="depth 2: value of a 2-hop link relative to a direct one [0.5]")
    p.add_argument("--max-shell", type=int, default=300,
                   help="depth 2: how many first-shell partners to expand [300]")
    p.add_argument("--cover", type=int, default=3, help="number of different targets to combine [3]")
    p.add_argument("--exact", action="store_true",
                   help="test all combinations of --cover targets (among --candidates) instead of greedy")
    p.add_argument("--candidates", type=int, default=60, help="--exact: candidate pool size [60]")
    p.add_argument("--mechanism", action="store_true",
                   help="annotate how proteins act on each other (activation/inhibition, "
                        "phosphorylation, ...) from OmniPath; human, mouse, rat")
    p.add_argument("--mech-targets", type=int, default=20,
                   help="how many top targets get path analysis and a view in network.html [20]")
    p.add_argument("--include-predicted", action="store_true",
                   help="--mechanism: also keep text-mined / predicted-only OmniPath records")
    p.add_argument("--omnipath-license", choices=["academic", "commercial"], default="academic",
                   help="OmniPath records allowed for this use [academic]")
    p.add_argument("--refresh-omnipath", action="store_true",
                   help="re-download the cached OmniPath tables (auto after 30 days)")
    p.add_argument("--structures", action="store_true",
                   help="look up PDB entries and bound ligands for the top targets and the "
                        "best pairs (RCSB + UniProt); adds structures.csv, ligands.csv, "
                        "pair_structures.csv")
    p.add_argument("--struct-targets", type=int, default=10,
                   help="--structures: how many top targets to look up [10]")
    p.add_argument("--pairs", type=int, default=3,
                   help="how many best target pairs to report (always printed and saved "
                        "to top_pairs.csv) [3]")
    p.add_argument("--ternary", action="store_true",
                   help="--structures: for each pair with a ligand on both sides, measure "
                        "ligand burial and exit vectors and scan clash-free ternary poses "
                        "(geometric triage for PROTACs - read the README caveats)")
    p.add_argument("--ternary-vectors", type=int, default=3,
                   help="--ternary: exit vectors per ligand to try [3]")
    p.add_argument("--ternary-poses", type=int, default=5,
                   help="--ternary: distinct arrangements to write as PDB files per pair [5]")
    p.add_argument("--refresh-structures", action="store_true",
                   help="re-download cached PDB metadata (auto after 30 days)")
    p.add_argument("--max-external", type=int, default=25, help="max outside proteins drawn [25]")
    p.add_argument("--top", type=int, default=10, help="rows to print [10]")
    p.add_argument("--out", default="results", help="output folder [results]")
    return p.parse_args()


def read_inputs(args):
    names = list(args.proteins)
    if args.file:
        for line in Path(args.file).read_text().splitlines():
            names += line.split("#")[0].replace(",", " ").split()
    return list(dict.fromkeys(n.strip() for n in names if n.strip()))


def fetch(args, use_string, use_biogrid, key, string_ids, genes, include_interactors):
    edges = []
    if use_string and string_ids:
        edges += string_api.interaction_partners(
            string_ids, args.species, args.score, args.network_type, args.limit, args.min_escore)
    if use_biogrid and genes:
        known = set(genes) | {n for e in edges for n in (e["a"], e["b"])}
        edges += biogrid_api.interactions(
            genes, key, args.biogrid_taxid or args.species,
            physical_only=not args.all_biogrid, min_publications=args.min_pubs,
            include_interactors=include_interactors, known_names=known)
    return edges


def run_ternary(args, pair_df, out):
    """Ligand burial, exit vectors and clash-free pose scan for each druggable pair."""
    from ppi import geometry

    usable = pair_df[pair_df.both_have_ligand]
    if not len(usable):
        print("\n--ternary: no pair has a ligand-bound structure on both sides, skipped")
        return []
    print(f"\nTernary geometry for {len(usable)} pair(s) "
          f"(rigid triage, not a ternary-complex prediction):")
    rows, vecs, poses = [], [], []
    for r in usable.itertuples():
        tag = f"{r.protein_a}_{r.protein_b}"
        try:
            summary, v, p = geometry.analyse_pair(
                r.ligand_pdb_a, r.ligand_pdb_b, n_vectors=args.ternary_vectors,
                out_dir=out / "poses", tag=tag, max_poses=args.ternary_poses)
        except Exception as exc:                                       # noqa: BLE001
            print(f"  {tag}: skipped ({exc})")
            continue
        summary = {"protein_a": r.protein_a, "protein_b": r.protein_b, **summary}
        rows.append(summary)
        v.insert(0, "pair", tag)
        vecs.append(v)
        p.insert(0, "pair", tag)
        poses.append(p.head(200))
        print(f"  {tag}: {r.ligand_pdb_a}/{summary['ligand_a']} buried "
              f"{summary['buried_fraction_a']:.0%}, {r.ligand_pdb_b}/{summary['ligand_b']} "
              f"buried {summary['buried_fraction_b']:.0%} -> "
              f"{summary['clash_free_pct']}% of {summary['n_poses']} poses clash-free, "
              f"{summary['distinct_arrangements']} distinct arrangements")
    if not rows:
        return []
    pd.DataFrame(rows).to_csv(out / "ternary_summary.csv", index=False)
    pd.concat(vecs, ignore_index=True).to_csv(out / "exit_vectors.csv", index=False)
    pd.concat(poses, ignore_index=True).to_csv(out / "ternary_poses.csv", index=False)
    print(f"  pose PDB files in {out}/poses/ - open them in PyMOL to inspect")
    return ["ternary_summary.csv", "exit_vectors.csv", "ternary_poses.csv", "poses/"]


def main():
    args = parse()
    names = read_inputs(args)
    if len(names) < 2:
        sys.exit("Give at least 2 proteins.")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    use_string = args.source in ("string", "both")
    use_biogrid = args.source in ("biogrid", "both")
    key = biogrid_api.get_key(args.biogrid_key) if use_biogrid else None

    # 1. Resolve names (STRING gives canonical gene names that BioGRID also uses)
    inputs, string_ids = list(names), []
    if use_string:
        print(f"Resolving {len(names)} names in STRING (species {args.species})...")
        res = string_api.resolve(names, args.species)
        inputs = []
        for n in names:
            if n in res:
                pref = res[n]["name"]
                if pref != n:
                    print(f"  {n} -> {pref}")
                if pref not in inputs:
                    inputs.append(pref)
                    string_ids.append(res[n]["string_id"])
            else:
                print(f"  not found in STRING: {n}")
                if use_biogrid:
                    inputs.append(n)
    canon = {n.upper(): n for n in inputs}
    inp = set(inputs)
    excluded = set()
    if args.exclude_hubs is not None:
        excluded = network.PROMISCUOUS | {n.upper() for n in args.exclude_hubs}
        kept = sorted(n for n in inputs if n.upper() in excluded)
        if kept:
            print(f"  not excluded because they are in your list: {', '.join(kept)}")

    # 2. Direct partners of the inputs
    print("Fetching interaction partners...")
    edges = fetch(args, use_string, use_biogrid, key, string_ids, inputs,
                  include_interactors=not args.only_within)
    edges = network.dedupe_edges(edges, canon)
    if excluded:
        edges = network.drop_proteins(edges, excluded, keep=inp)
    if args.only_within:
        edges = [e for e in edges if e["a"] in inp and e["b"] in inp]
    if not edges:
        sys.exit("No interactions found - lower --score / --min-escore or check names and species.")
    counts = network.partner_counts(edges, inputs)

    # 3. Depth 2: links among the first shell, so outside proteins can bridge to inputs
    if args.depth == 2 and not args.only_within:
        G1 = network.build_graph(edges, inputs)
        r1 = network.rank_targets(G1, inputs)
        shell = r1[(~r1.in_input) & (r1.n_input_partners >= 1)].protein.head(args.max_shell).tolist()
        print(f"Depth 2: fetching links among {len(shell)} first-shell proteins...")
        id_of = {}
        for e in edges:
            if e["source"] == "STRING":
                id_of[e["a"]], id_of[e["b"]] = e["id_a"], e["id_b"]
        new = fetch(args, use_string, use_biogrid, key,
                    [id_of[n] for n in shell if n in id_of], inputs + shell,
                    include_interactors=False)
        keep = inp | set(shell)
        new = [e for e in network.dedupe_edges(new, canon) if e["a"] in keep and e["b"] in keep]
        new = network.drop_proteins(new, excluded, keep=inp)
        edges = network.dedupe_edges(edges + new)

    # 4. Analyse
    pd.DataFrame(edges).to_csv(out / "edges.csv", index=False)
    counts.to_csv(out / "input_partner_counts.csv", index=False)
    G = network.build_graph(edges, inputs)
    table = network.reach_table(G, inputs, args.depth)
    ranking = network.rank_targets(G, inputs, table, args.depth, args.hop_weight)
    cover_fn = network.exact_cover if args.exact else network.greedy_cover
    kw = {"candidates": args.candidates} if args.exact else {}
    picks, left = cover_fn(table, inputs, args.cover, args.hop_weight, **kw)
    pd.DataFrame(picks).to_csv(out / "cover.csv", index=False)

    # 5. Mechanisms: how the top targets act on the inputs (directly or via a bridge)
    index, saved = {}, ["target_ranking.csv", "cover.csv", "input_partner_counts.csv", "edges.csv"]
    if args.mechanism:
        if args.species in mechanisms.SPECIES:
            print("Annotating mechanisms (OmniPath)...")
            inter, enz = mechanisms.load_tables(args.species, args.omnipath_license,
                                                args.refresh_omnipath)
            index = mechanisms.build_index(inter, enz, G.nodes, args.include_predicted)
            mechanisms.index_table(index, G).to_csv(out / "mechanisms.csv", index=False)
            saved.append("mechanisms.csv")
        else:
            print(f"--mechanism: OmniPath has no data for species {args.species}, skipped")
    targets = list(dict.fromkeys(ranking.protein.head(args.mech_targets).tolist()
                                 + [p["target"] for p in picks]))
    paths_df, effects, views = mechanisms.analyse_targets(G, index, table, targets)
    if index:
        paths_df.to_csv(out / "paths.csv", index=False)
        effects.to_csv(out / "target_effects.csv", index=False)
        saved += ["paths.csv", "target_effects.csv"]
        ranking = ranking.merge(effects[[c for c in effects if c == "protein" or c.startswith("n_")]],
                                on="protein", how="left")
    ranking.to_csv(out / "target_ranking.csv", index=False)

    # 6. Structures: what is known in the PDB, and can a pair be linked without clashing
    pairs = network.top_pairs(table, inputs, args.pairs, args.hop_weight)
    pd.DataFrame(pairs).to_csv(out / "top_pairs.csv", index=False)
    saved.append("top_pairs.csv")
    if args.structures:
        want = list(dict.fromkeys(ranking.protein.head(args.struct_targets).tolist()
                                  + [p[k] for p in pairs for k in ("protein_a", "protein_b")]
                                  + inputs))
        print(f"Looking up PDB structures for {len(want)} proteins...")
        sdf, accs, meta, per_gene = structures.collect(want, args.species,
                                                       args.refresh_structures)
        if len(sdf):
            sdf.to_csv(out / "structures.csv", index=False)
            structures.ligand_table(meta, {a: g for g, a in accs.items()}).to_csv(
                out / "ligands.csv", index=False)
            pdf = structures.pair_table([(p["protein_a"], p["protein_b"]) for p in pairs],
                                        accs, per_gene, meta)
            pdf.to_csv(out / "pair_structures.csv", index=False)
            saved += ["structures.csv", "ligands.csv", "pair_structures.csv"]
            print("\nStructures available:")
            print(sdf.head(args.top)[["protein", "uniprot", "n_structures", "n_with_ligand",
                                      "best_resolution", "best_ligand_pdb",
                                      "n_distinct_ligands"]].to_string(index=False))
            print(f"\nBest {len(pairs)} target pairs:")
            print(pdf[["protein_a", "protein_b", "ligand_pdb_a", "ligand_pdb_b",
                       "n_cocomplex", "cocomplex_pdbs", "both_have_ligand"]]
                  .to_string(index=False))
            if args.ternary:
                saved += run_ternary(args, pdf, out)
        else:
            print("  no structures found for these proteins")

    # 7. Draw
    best = ranking.iloc[0].protein
    paths = network.two_hop_paths(G, best, table) if args.depth == 2 else []
    extra = [p["target"] for p in picks] + [x for x, _ in paths]
    H = plot.focus_subgraph(G, ranking, args.max_external, extra=extra)
    plot.draw_png(H, ranking, out / "network.png", best, paths)
    interactive.draw_html(G, ranking, out / "network.html", set(H.nodes), views, index, best,
                          title=f"{len(inputs)} proteins: {', '.join(inputs[:8])}"
                                + (" ..." if len(inputs) > 8 else ""))
    saved += ["network.png", "network.html"]

    # 8. Report
    cols = ["protein", "in_input", "n_input_partners", "coverage_pct", "input_partners"]
    if args.depth == 2:
        cols[3:3] = ["n_two_hop", "score"]
    with pd.option_context("display.max_colwidth", 70, "display.width", 220):
        print(f"\n{len(inputs)} input proteins | {G.number_of_nodes()} proteins and "
              f"{G.number_of_edges()} interactions | depth {args.depth}")
        print("\nPartners per input protein:")
        print(counts.to_string(index=False))
        print("\nBest single targets:")
        print(ranking.head(args.top)[cols].to_string(index=False))
        if pairs:
            print(f"\nBest {len(pairs)} target pairs (each protein used once):")
            print(pd.DataFrame(pairs).to_string(index=False))
        if picks:
            how = f"exact, best of top {args.candidates}" if args.exact else "greedy"
            print(f"\nBest combination of {len(picks)} different targets ({how}):")
            print(pd.DataFrame(picks).to_string(index=False))
        if left:
            print("Not reached by these targets:", ", ".join(left))
        if index:
            print("\nHow the top targets act on your list (OmniPath; details in paths.csv):")
            cols = ["activates", "inhibits", "mixed", "sign_unknown", "no_mechanism_data"]
            e = effects.head(args.top)[["protein"] + [f"n_{c}" for c in cols]]
            print(e.rename(columns=lambda c: c[2:] if c.startswith("n_") else c)
                  .to_string(index=False))
    print(f"\nSaved to {out}/: {', '.join(saved)}")


if __name__ == "__main__":
    main()
