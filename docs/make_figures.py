#!/usr/bin/env python3
"""Regenerate the presentation figures in docs/.

    python docs/make_figures.py

Writes pipeline.png, ternary_scan.png and benchmark_5t35.png at 200 dpi on a light
surface, sized for a 16:9 slide. Colours come from a palette checked for
colour-blind separation; every element is directly labelled, so nothing depends on
colour alone.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle, Wedge, Polygon  # noqa: E402

OUT = Path(__file__).resolve().parent
SURFACE = "#fcfcfb"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8a85"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
FONT = {"family": "DejaVu Sans"}
plt.rcParams.update({"font.family": FONT["family"], "savefig.facecolor": SURFACE})


def _box(ax, x, y, w, h, text, face, edge=None, fs=11, weight="bold", tc="white"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.02",
                                linewidth=1.2 if edge else 0, edgecolor=edge or "none",
                                facecolor=face, zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", color=tc,
            fontsize=fs, fontweight=weight, zorder=3, linespacing=1.35)


def _arrow(ax, p0, p1, color=INK3, lw=1.6, style="-|>", rad=0.0):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=13,
                                 linewidth=lw, color=color, zorder=1,
                                 connectionstyle=f"arc3,rad={rad}"))


# ---------------------------------------------------------------- pipeline

def pipeline():
    fig, ax = plt.subplots(figsize=(13.33, 7.5))
    fig.patch.set_facecolor(SURFACE)
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

    ax.text(0.15, 9.6, "ppi-hub-finder", fontsize=19, fontweight="bold", color=INK)
    ax.text(0.15, 9.18, "From a list of proteins to a degrader hypothesis",
            fontsize=12.5, color=INK2)

    ax.text(1.55, 8.55, "DATA SOURCES", fontsize=8.5, color=INK3,
            fontweight="bold", ha="center")
    ax.text(5.0, 8.55, "PIPELINE", fontsize=8.5, color=INK3, fontweight="bold", ha="center")
    ax.text(8.6, 8.55, "WHAT YOU GET", fontsize=8.5, color=INK3,
            fontweight="bold", ha="center")

    stages = [
        ("1. Resolve names\nSTRING canonical symbols", "STRING", "", 7.55),
        ("2. Fetch interactions\ndirect partners, optional 2nd step",
         "STRING\nBioGRID", "edges.csv", 6.15),
        ("3. Rank targets and pairs\nwho reaches most of your list",
         "", "target_ranking.csv\ntop_pairs.csv\nnetwork.html", 4.75),
        ("4. Mechanism\nactivates / inhibits, which residue",
         "OmniPath", "paths.csv\ntarget_effects.csv", 3.35),
        ("5. Structures and ligands\nper protein, per pair, per E3",
         "RCSB PDB\nUniProt", "structures.csv\ne3_pairs.csv", 1.95),
        ("6. Ternary geometry\nburial, exit vectors, clash-free poses",
         "built-in\nE3 table", "ternary_summary.csv\npose_*.pdb", 0.55),
    ]
    bw, bh, bx = 3.5, 1.0, 3.25
    for i, (label, src, out, y) in enumerate(stages):
        head, sub = label.split("\n")
        _box(ax, bx, y, bw, bh, "", BLUE if i < 3 else ORANGE)
        ax.text(bx + 0.22, y + bh * 0.63, head, fontsize=11.5, fontweight="bold",
                color="white", va="center")
        ax.text(bx + 0.22, y + bh * 0.26, sub, fontsize=9.2, color="#f2f2ef", va="center")
        if i:
            _arrow(ax, (bx + bw / 2, y + bh + 0.38), (bx + bw / 2, y + bh + 0.04))
        if src:
            _box(ax, 0.45, y + 0.12, 2.2, bh - 0.24, src, "none", edge=INK3,
                 fs=9.5, weight="normal", tc=INK2)
            _arrow(ax, (2.68, y + bh / 2), (bx - 0.03, y + bh / 2), color=INK3, lw=1.2)
        if out:
            _arrow(ax, (bx + bw + 0.03, y + bh / 2), (bx + bw + 0.42, y + bh / 2),
                   color=AQUA, lw=1.4)
            ax.text(bx + bw + 0.55, y + bh / 2, out, fontsize=9.2, color=INK2,
                    va="center", linespacing=1.4)

    ax.text(bx + 0.22, 0.12, "Steps 1-3 always run. 4, 5 and 6 are opt-in flags "
            "(--mechanism, --structures, --e3 --ternary).",
            fontsize=9, color=INK3, style="italic")
    fig.savefig(OUT / "pipeline.png", dpi=200, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


# ---------------------------------------------------------- ternary scan

def ternary_scan():
    fig, ax = plt.subplots(figsize=(13.33, 6.4))
    fig.patch.set_facecolor(SURFACE)
    ax.set_xlim(0, 10); ax.set_ylim(0, 6.4); ax.axis("off")

    ax.text(0.2, 6.05, "How the ternary scan works", fontsize=18,
            fontweight="bold", color=INK)
    ax.text(0.2, 5.68, "Protein A is held still. Protein B is placed so the two ligand "
            "exit vectors face each other, then four parameters are swept.",
            fontsize=11.5, color=INK2)

    ya = 3.45
    # --- protein A (fixed)
    ax.add_patch(Circle((1.55, ya), 1.05, facecolor=BLUE, alpha=0.16,
                        edgecolor=BLUE, linewidth=1.6, zorder=1))
    ax.text(1.55, ya + 1.32, "Protein A  (fixed)", ha="center", fontsize=11.5,
            fontweight="bold", color=BLUE)
    ax.add_patch(Circle((2.28, ya), 0.29, facecolor=BLUE, zorder=3))
    ax.text(2.28, ya, "L", ha="center", va="center", color="white",
            fontsize=10, fontweight="bold", zorder=4)

    # --- exit vector
    _arrow(ax, (2.60, ya), (3.85, ya), color=ORANGE, lw=2.4)
    ax.text(3.22, ya + 0.30, "exit vector", ha="center", fontsize=10.5,
            fontweight="bold", color=ORANGE)
    ax.text(2.28, ya - 0.52, "bound ligand", ha="center", fontsize=9, color=INK2)
    ax.text(2.28, ya - 0.84, "its most solvent-exposed atom,\npointing away from the pocket",
            ha="center", fontsize=8.6, color=INK3, va="top", linespacing=1.35)

    # --- tilt fan (drawn faint, labelled once)
    for ang, sign in ((0.40, 1), (0.40, -1)):
        _arrow(ax, (2.60, ya), (2.60 + 1.25, ya + sign * ang), color=AQUA, lw=1.1)
    ax.text(3.95, ya + 0.62, "tilt", fontsize=9.5, fontweight="bold", color=AQUA)

    # --- protein B (moved)
    ax.add_patch(Circle((6.05, ya), 1.12, facecolor=ORANGE, alpha=0.16,
                        edgecolor=ORANGE, linewidth=1.6, zorder=1))
    ax.text(6.05, ya + 1.40, "Protein B or E3  (moved)", ha="center", fontsize=11.5,
            fontweight="bold", color=ORANGE)
    ax.add_patch(Circle((5.12, ya), 0.29, facecolor=ORANGE, zorder=3))
    ax.text(5.12, ya, "L", ha="center", va="center", color="white",
            fontsize=10, fontweight="bold", zorder=4)
    # spin indicator around B
    ax.add_patch(FancyArrowPatch((6.95, ya + 0.72), (6.95, ya - 0.72),
                                 arrowstyle="-|>", mutation_scale=12, linewidth=1.5,
                                 color=INK3, connectionstyle="arc3,rad=-0.55", zorder=3))
    ax.text(7.32, ya, "spin", fontsize=9.5, fontweight="bold", color=INK3, va="center")

    # --- linker span
    ax.annotate("", xy=(4.83, ya - 1.45), xytext=(2.57, ya - 1.45),
                arrowprops=dict(arrowstyle="<|-|>", color=INK, lw=1.4))
    ax.text(3.70, ya - 1.78, "linker span", ha="center", fontsize=10,
            fontweight="bold", color=INK)

    # --- the four swept parameters, listed rather than crowded onto the drawing
    _box(ax, 7.95, ya - 1.62, 1.95, 3.02, "", "#f1f1ee")
    ax.text(8.12, ya + 1.18, "SWEPT", fontsize=8.5, fontweight="bold", color=INK3)
    rows = [("1", "linker span", "5-20 A", INK),
            ("2", "tilt", "0 / 20 / 40 deg", AQUA),
            ("3", "azimuth", "6 directions", INK3),
            ("4", "spin", "24 angles", INK3)]
    for i, (n, name, detail, col) in enumerate(rows):
        yy = ya + 0.80 - i * 0.70
        ax.add_patch(Circle((8.28, yy), 0.145, facecolor=col, zorder=3))
        ax.text(8.28, yy, n, ha="center", va="center", color="white",
                fontsize=8.5, fontweight="bold", zorder=4)
        ax.text(8.52, yy + 0.09, name, fontsize=10, fontweight="bold", color=INK)
        ax.text(8.52, yy - 0.19, detail, fontsize=8.8, color=INK2)
    ax.text(8.12, ya - 1.92, "about 17 000 poses per vector pair",
            fontsize=8.8, color=INK3)

    # --- scoring strip
    _box(ax, 0.2, 0.18, 9.7, 1.12, "", "#f1f1ee")
    ax.text(0.48, 1.02, "Every pose is scored", fontsize=11, fontweight="bold", color=INK)
    ax.text(0.48, 0.70, "clash = heavy atoms closer than 2.5 A     contact = within 4.5 A"
            "     near-identical poses are clustered away",
            fontsize=9.8, color=INK2)
    ax.text(0.48, 0.38, "Reported: what fraction of poses are clash-free, and a handful of "
            "genuinely different arrangements as PDB files.", fontsize=9.8, color=INK2)
    fig.savefig(OUT / "ternary_scan.png", dpi=200, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


# ------------------------------------------------------------- benchmark

def benchmark():
    """What the 5T35 check showed: sampling reaches the native pose, scoring misses it."""
    fig, ax = plt.subplots(figsize=(11, 5.4))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    fig.subplots_adjust(left=0.30, right=0.96, top=0.80, bottom=0.26)

    labels = ["Closest pose the scan\nproduced  (sampling)",
              "Pose our score ranked\nfirst  (scoring)",
              "Median clash-free pose"]
    values = [5.4, 30.6, 30.8]
    colors = [AQUA, ORANGE, ORANGE]
    ax.barh(range(len(values)), values, height=0.5, color=colors, zorder=3)
    for i, v in enumerate(values):
        ax.text(v + 0.6, i, f"{v} A", va="center", fontsize=12.5,
                fontweight="bold", color=INK)

    ax.axvline(5.0, color=INK3, linestyle=(0, (4, 3)), linewidth=1.3, zorder=2)
    ax.text(5.35, -0.46, "about the size of a native-like arrangement",
            fontsize=9, color=INK3, va="center")

    ax.set_yticks(range(len(values)))
    ax.set_yticklabels(labels, fontsize=11, color=INK, linespacing=1.4)
    ax.set_ylim(2.6, -0.7)
    ax.set_xlabel("CA RMSD to the crystal structure of the real complex (A)  -  lower is better",
                  fontsize=10.5, color=INK2, labelpad=9)
    ax.set_xlim(0, 36)
    ax.tick_params(axis="x", colors=INK2, labelsize=10)
    ax.tick_params(axis="y", length=0)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(INK3)
    ax.grid(axis="x", color="#e6e6e2", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)

    fig.text(0.02, 0.93, "Validation against BRD4 : MZ1 : VHL  (PDB 5T35)",
             fontsize=15.5, fontweight="bold", color=INK)
    fig.text(0.02, 0.875, "Rebuilt from the two binary structures the scan would "
             "normally start from, then compared with the real complex.",
             fontsize=10.5, color=INK2)
    fig.text(0.02, 0.055, "The native arrangement IS among the poses generated - but "
             "ranking by clashes and contacts does not find it.",
             fontsize=10.2, color=INK, fontweight="bold")
    fig.text(0.02, 0.012, "So read the output as triage (can these two be linked at "
             "all?), not as a predicted ternary complex.", fontsize=10.2, color=INK2)
    fig.savefig(OUT / "benchmark_5t35.png", dpi=200, facecolor=SURFACE)
    plt.close(fig)


if __name__ == "__main__":
    pipeline()
    ternary_scan()
    benchmark()
    print("wrote", ", ".join(sorted(p.name for p in OUT.glob("*.png"))))
