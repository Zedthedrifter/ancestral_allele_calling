#!/usr/bin/env python3
"""
Plot genome-wide density of ancestral allele disagreements.
Usage: python3 plot_disagreements.py disagreements.tsv out_prefix
"""
import sys
import numpy as np
import matplotlib.pyplot as plt

infile = sys.argv[1]
prefix = sys.argv[2]

# ---- Load data ----
chrom, pos = [], []
with open(infile) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        chrom.append(f[0])
        pos.append(int(f[1]))

chrom = np.array(chrom)
pos   = np.array(pos)

# ---- Chromosome order (natural sort for 1..12) ----
chroms = sorted(set(chrom), key=lambda x: int(''.join(filter(str.isdigit, x))))

# ---- Bin size ----
BIN = 100_000   # 100 kb windows

fig, axes = plt.subplots(len(chroms), 1,
                         figsize=(12, 1.6 * len(chroms)),
                         sharex=False)

for ax, c in zip(axes, chroms):
    m = chrom == c
    p = pos[m]
    maxpos = p.max() if len(p) else 1
    bins = np.arange(0, maxpos + BIN, BIN)
    counts, edges = np.histogram(p, bins=bins)
    centers = edges[:-1] + BIN / 2

    ax.bar(centers, counts, width=BIN * 0.9,
           color="#3b6ea5", edgecolor="none")

    # ---- Y-axis: show integer ticks and label ----
    ymax = counts.max() if len(counts) else 1
    ax.set_ylim(0, ymax * 1.15)
    ax.set_yticks(np.arange(0, ymax + 1, max(1, ymax // 5)))
    ax.tick_params(axis="y", labelsize=8)
    ax.set_ylabel('count', rotation=0, ha="right", va="center", labelpad=12)

    # ---- Chromosome label ----
    ax.text(1.005, 0.5, c, transform=ax.transAxes,
            ha="left", va="center", fontsize=9)

    ax.set_xlim(0, maxpos)

axes[-1].set_xlabel(f"Position (bins of {BIN:,} bp)")
fig.suptitle("Ancestral allele disagreements across the genome", y=1.0)
fig.tight_layout(rect=[0, 0, 0.97, 1.0])
fig.savefig(f"{prefix}_density.png", dpi=200)
print(f"Saved {prefix}_density.png")

# ---- Also write the binned table for downstream use ----
with open(f"{prefix}_density.tsv", "w") as fh:
    fh.write("chrom\twindow_start\twindow_end\tn_disagreements\n")
    for c in chroms:
        m = chrom == c
        p = pos[m]
        maxpos = p.max() if len(p) else 1
        bins = np.arange(0, maxpos + BIN, BIN)
        counts, edges = np.histogram(p, bins=bins)
        for i, ct in enumerate(counts):
            fh.write(f"{c}\t{int(edges[i])}\t{int(edges[i+1])}\t{ct}\n")
print(f"Saved {prefix}_density.tsv")