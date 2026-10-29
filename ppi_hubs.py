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

from ppi import biogrid_api, network, plot, string_api


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
    p.add_argument("--min-pubs", type=int, default=1, help="BioGRID min publications per pair [1]")
    p.add_argument("--all-biogrid", action="store_true", help="also keep BioGRID genetic interactions")
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
        edges += biogrid_api.interactions(
            genes, key, args.biogrid_taxid or args.species,
            physical_only=not args.all_biogrid, min_publications=args.min_pubs,
            include_interactors=include_interactors)
    return edges


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

    # 2. Direct partners of the inputs
    print("Fetching interaction partners...")
    edges = fetch(args, use_string, use_biogrid, key, string_ids, inputs,
                  include_interactors=not args.only_within)
    edges = network.dedupe_edges(edges, canon)
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
                    [id_of[n] for n in shell if n in id_of], shell, include_interactors=True)
        keep = inp | set(shell)
        new = [e for e in network.dedupe_edges(new, canon) if e["a"] in keep and e["b"] in keep]
        edges = network.dedupe_edges(edges + new)

    # 4. Analyse
    pd.DataFrame(edges).to_csv(out / "edges.csv", index=False)
    counts.to_csv(out / "input_partner_counts.csv", index=False)
    G = network.build_graph(edges, inputs)
    table = network.reach_table(G, inputs, args.depth)
    ranking = network.rank_targets(G, inputs, table, args.depth, args.hop_weight)
    ranking.to_csv(out / "target_ranking.csv", index=False)
    cover_fn = network.exact_cover if args.exact else network.greedy_cover
    kw = {"candidates": args.candidates} if args.exact else {}
    picks, left = cover_fn(table, inputs, args.cover, args.hop_weight, **kw)
    pd.DataFrame(picks).to_csv(out / "cover.csv", index=False)

    # 5. Draw
    best = ranking.iloc[0].protein
    paths = network.two_hop_paths(G, best, table) if args.depth == 2 else []
    extra = [p["target"] for p in picks] + [x for x, _ in paths]
    H = plot.focus_subgraph(G, ranking, args.max_external, extra=extra)
    plot.draw_png(H, ranking, out / "network.png", best, paths)
    plot.draw_html(H, ranking, out / "network.html", best, paths)

    # 6. Report
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
        if picks:
            how = f"exact, best of top {args.candidates}" if args.exact else "greedy"
            print(f"\nBest combination of {len(picks)} different targets ({how}):")
            print(pd.DataFrame(picks).to_string(index=False))
        if left:
            print("Not reached by these targets:", ", ".join(left))
    print(f"\nSaved to {out}/: target_ranking.csv, cover.csv, input_partner_counts.csv, "
          "edges.csv, network.png, network.html")


if __name__ == "__main__":
    main()
