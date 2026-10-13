#!/usr/bin/env python3
"""
Plot disagreement positions as points on chromosome ideograms.
Usage: python3 plot_disagreements_points.py disagreements.tsv out_prefix
"""
import sys
import numpy as np
import matplotlib.pyplot as plt

infile = sys.argv[1]
prefix = sys.argv[2]

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
chroms = sorted(set(chrom), key=lambda x: int(''.join(filter(str.isdigit, x))))

# ---- Get chromosome lengths from the data (or from a reference) ----
# Here we use the max observed position per chromosome as a proxy.
lengths = {c: pos[chrom == c].max() for c in chroms}

fig, ax = plt.subplots(figsize=(10, 6))

for i, c in enumerate(chroms):
    y = len(chroms) - i
    L = lengths[c]
    # Chromosome backbone
    ax.plot([0, L], [y, y], color="lightgrey", lw=6, solid_capstyle="butt")
    # Disagreement points
    p = pos[chrom == c]
    ax.scatter(p, np.full_like(p, y), s=4, color="#a5523b", alpha=0.6)

ax.set_yticks([len(chroms) - i for i in range(len(chroms))])
ax.set_yticklabels(chroms)
ax.set_xlabel("Position (bp)")
ax.set_title("Positions of ancestral allele disagreements")
ax.set_ylim(0.5, len(chroms) + 0.5)
fig.tight_layout()
fig.savefig(f"{prefix}_points.png", dpi=200)
print(f"Saved {prefix}_points.png")
