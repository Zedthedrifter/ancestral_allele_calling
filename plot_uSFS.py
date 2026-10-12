#!/usr/bin/env python3
"""
Compute and plot the unfolded site frequency spectrum (uSFS) from an
est-sfs ancestral allele TSV.

Expected input columns (tab-delimited):
  1 chrom
  2 pos
  3 ancestral allele  (A/C/G/T or NA)
  4 major allele
  5 P-major-ancestral
  6 P(A)
  7 P(C)
  8 P(G)
  9 P(T)
 10 focal counts (A,C,G,T)
 11 OG1 state     (A,C,G,T)
 12 OG2 state     (A,C,G,T)

Usage:
  python3 plot_uSFS.py input.tsv output_prefix
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Arguments
# ----------------------------------------------------------------------
infile  = sys.argv[1]
prefix  = sys.argv[2]
outdir  = sys.argv[3]

# ----------------------------------------------------------------------
# Read the file and compute derived allele counts
# ----------------------------------------------------------------------
NUC_TO_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}

derived_counts = []   # number of derived alleles at each site
n_total        = []   # total number of alleles at each site (n)
n_skipped_aa   = 0
n_skipped_zero = 0

with open(infile) as fh:
    for line in fh:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 10:
            continue

        aa = fields[2]
        if aa not in NUC_TO_IDX:
            n_skipped_aa += 1
            continue

        focal = fields[9].split(",")
        if len(focal) != 4:
            continue

        counts = [int(x) for x in focal]
        n = sum(counts)
        if n == 0:
            n_skipped_zero += 1
            continue

        # Derived = non-ancestral
        anc_idx = NUC_TO_IDX[aa]
        derived = n - counts[anc_idx]

        derived_counts.append(derived)
        n_total.append(n)

derived_counts = np.array(derived_counts)
n_total        = np.array(n_total)

print(f"Sites used: {len(derived_counts)}")
print(f"Skipped (NA ancestral): {n_skipped_aa}")
print(f"Skipped (zero counts):  {n_skipped_zero}")

# ----------------------------------------------------------------------
# Build the uSFS
# ----------------------------------------------------------------------
# For est-sfs, n is the constant --focalN used (e.g. 100), but the raw
# allele counts may differ per site if the focal file was not down-sampled.
# We take the mode of n_total as the canonical n.
n_mode = int(np.bincount(n_total).argmax())
print(f"Modal number of alleles per site: {n_mode}")

# Filter to sites with exactly n_mode alleles
mask = (n_total == n_mode)
dc   = derived_counts[mask]

# uSFS: number of sites with k derived alleles, k = 1 .. n_mode - 1
usfs = np.zeros(n_mode, dtype=int)
for k in dc:
    if 1 <= k <= n_mode - 1:
        usfs[k] += 1

# ----------------------------------------------------------------------
# Save uSFS as a table
# ----------------------------------------------------------------------
out_table = f"{outdir}/{prefix}_uSFS.tsv"
with open(out_table, "w") as fh:
    fh.write("derived_count\tn_sites\n")
    for k in range(1, n_mode):
        fh.write(f"{k}\t{usfs[k]}\n")
print(f"uSFS table written to {out_table}")

# ----------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 4.5))
k = np.arange(1, n_mode)
ax.bar(k, usfs[1:n_mode], width=0.9, color="#3b6ea5", edgecolor="none")
ax.set_xlabel("Derived allele count")
ax.set_ylabel("Number of sites")
ax.set_title(f"Unfolded SFS (n = {n_mode})")
ax.set_xlim(0, n_mode)
fig.tight_layout()

out_fig = f"{outdir}/{prefix}_uSFS.png"
fig.savefig(out_fig, dpi=200)
print(f"uSFS plot written to {out_fig}")

# Optional: also save a log-scaled version (useful for the tails)
fig2, ax2 = plt.subplots(figsize=(7, 4.5))
ax2.bar(k, usfs[1:n_mode], width=0.9, color="#a5523b", edgecolor="none")
ax2.set_yscale("log")
ax2.set_xlabel("Derived allele count")
ax2.set_ylabel("Number of sites (log scale)")
ax2.set_title(f"Unfolded SFS, log scale (n = {n_mode})")
ax2.set_xlim(0, n_mode)
fig2.tight_layout()

out_fig_log = f"{outdir}/{prefix}_uSFS_log.png"
fig2.savefig(out_fig_log, dpi=200)
print(f"Log-scale uSFS plot written to {out_fig_log}")
