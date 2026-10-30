"""Static (PNG) interaction map. The interactive map is in interactive.py."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

INPUT_C, HUB_C, BEST_C, EDGE_C = "#4C72B0", "#DD8452", "#C44E52", "#B0B0B0"


def focus_subgraph(G, ranking, max_external=25, min_score=2, extra=()):
    """Input proteins + the best outside proteins (score >= min_score) + any extra nodes."""
    inputs = [n for n, d in G.nodes(data=True) if d["is_input"]]
    ext = ranking[(~ranking.in_input) & (ranking.score >= min_score)]
    keep = set(inputs) | set(ext.protein.head(max_external)) | set(extra)
    return G.subgraph(keep).copy()


def _split_edges(H, best, paths):
    path_set = {frozenset(p) for p in paths}
    best_e = [e for e in H.edges if best in e]
    path_e = [e for e in H.edges if best not in e and frozenset(e) in path_set]
    other = [e for e in H.edges if best not in e and frozenset(e) not in path_set]
    return best_e, path_e, other


def draw_png(H, ranking, path, best=None, paths=()):
    info = ranking.set_index("protein")
    is_in = nx.get_node_attributes(H, "is_input")
    pos = nx.spring_layout(H, seed=42, k=2.0 / max(len(H), 1) ** 0.5, iterations=200)
    fig, ax = plt.subplots(figsize=(14, 11))

    best_e, path_e, other = _split_edges(H, best, paths)
    nx.draw_networkx_edges(H, pos, edgelist=other, edge_color=EDGE_C, ax=ax,
                           width=[0.6 + 1.6 * (len(H.edges[e]["sources"]) - 1) for e in other])
    nx.draw_networkx_edges(H, pos, edgelist=best_e, width=2.4, edge_color=BEST_C, ax=ax)
    if path_e:
        nx.draw_networkx_edges(H, pos, edgelist=path_e, width=1.8, edge_color=BEST_C,
                               style="dashed", ax=ax)

    for group, color, shape in ((True, INPUT_C, "o"), (False, HUB_C, "D")):
        nodes = [n for n in H if is_in[n] == group]
        if not nodes:
            continue
        nx.draw_networkx_nodes(
            H, pos, nodelist=nodes, node_color=color, node_shape=shape, ax=ax,
            node_size=[250 + 160 * float(info.at[n, "score"]) for n in nodes],
            edgecolors=[BEST_C if n == best else "white" for n in nodes],
            linewidths=[3.5 if n == best else 1.0 for n in nodes])
    nx.draw_networkx_labels(H, pos, font_size=8, ax=ax)

    handles = [
        Line2D([], [], marker="o", ls="", markerfacecolor=INPUT_C, markeredgecolor="w",
               markersize=11, label="Your input protein"),
        Line2D([], [], marker="D", ls="", markerfacecolor=HUB_C, markeredgecolor="w",
               markersize=10, label="Outside protein linking several inputs"),
        Line2D([], [], color=BEST_C, lw=2.4, label="Best target: direct links"),
        Line2D([], [], color=EDGE_C, lw=2.2, label="Thick grey = in both STRING and BioGRID"),
    ]
    if path_e:
        handles.insert(3, Line2D([], [], color=BEST_C, lw=1.8, ls="--",
                                 label="Best target: 2-hop paths (via intermediate)"))
    ax.legend(handles=handles, loc="lower left", fontsize=9, frameon=False)
    r = info.loc[best]
    two = f" + {int(r.n_two_hop)} via one intermediate" if r.n_two_hop else ""
    ax.set_title(f"Best target: {best} - {int(r.n_input_partners)} direct{two} "
                 f"of {sum(is_in.values())} input proteins (node size = score)", fontsize=12)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
