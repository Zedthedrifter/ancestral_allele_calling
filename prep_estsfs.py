#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Prepare an input file for Keightley & Jackson (2018) est-sfs.

For every biallelic SNP in a (optionally gzipped) VCF, writes one row:
  chrom  pos  ref  alt  MAJOR  MINOR  focal[A,C,G,T]

Column 7 (focal) is the est-sfs data file.
The leading columns are kept for bookkeeping / mapping est-sfs output
back to genomic coordinates and to your major/minor alleles.

Usage:
  python3 prep_estsfs.py \
      --infile in.vcf.gz \
      --focalSampleList focal.txt \
      --outfile focal.tsv \
      --focalN 100 --seed 42
"""
import gzip
import argparse
import re
import numpy as np

# ----------------------------------------------------------------------
# Arguments
# ----------------------------------------------------------------------
p = argparse.ArgumentParser(description="Generate est-sfs focal input from a VCF")
p.add_argument("--infile", required=True)
p.add_argument("--focalSampleList", required=True)
p.add_argument("--outfile", required=True)
# est-sfs needs a CONSTANT ingroup sample size. If set, every site is randomly
# down-sampled (without replacement) to this many alleles; sites with fewer
# called alleles are dropped. Strongly recommended for real data.
p.add_argument("--focalN", type=int, default=None,
               help="target number of focal ALLELES (constant n for est-sfs)")
p.add_argument("--seed", type=int, default=42,
               help="RNG seed for down-sampling")
p.add_argument("--keepMono", action="store_true",
               help="also write sites that are monomorphic in the focal group")
args = p.parse_args()

rng = np.random.default_rng(args.seed)

# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def smart_open(path):
    """Open plain or gzipped text transparently."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, "rt")

def read_list(path):
    with open(path) as fh:
        return fh.read().splitlines()

GT_SPLIT = re.compile(r"[/|]")

def count_ref_alt(calls):
    """
    Count called REF(0) and ALT(1) alleles across a list of GT fields.
    Robust to phasing, half-calls, haploid calls, and any missing form.
    """
    an = ac = 0
    for c in calls:
        gt = c.partition(":")[0]
        for a in GT_SPLIT.split(gt):
            if a == "0":
                an += 1
            elif a == "1":
                an += 1
                ac += 1
    return an, ac

NUCS = ["A", "C", "G", "T"]

# ----------------------------------------------------------------------
# Load sample list
# ----------------------------------------------------------------------
focal = read_list(args.focalSampleList)

# ----------------------------------------------------------------------
# Read VCF header once
# ----------------------------------------------------------------------
inVCF = smart_open(args.infile)
samples = None
for line in inVCF:
    if line.startswith("##"):
        continue
    if line.startswith("#CHROM"):
        samples = line.rstrip("\n").split("\t")[9:]
        break
if samples is None:
    raise SystemExit("No #CHROM header line found.")

sample_idx = {s: i for i, s in enumerate(samples)}
missing = [s for s in focal if s not in sample_idx]
if missing:
    raise SystemExit(f"Samples not in VCF: {missing[:10]} ...")

focal_cols = [sample_idx[s] for s in focal]

# ----------------------------------------------------------------------
# Main loop
# ----------------------------------------------------------------------
out = open(args.outfile, "w")
n_written = 0

for line in inVCF:
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    chrom, pos, ref, alt = f[0], f[1], f[3], f[4]

    if ref not in NUCS or alt not in NUCS:
        continue

    geno = f[9:]

    focal_calls = [geno[i] for i in focal_cols]
    AN, AC = count_ref_alt(focal_calls)
    if AN == 0:
        continue

    if args.focalN is not None:
        if AN < args.focalN:
            continue
        AC = int(rng.hypergeometric(AC, AN - AC, args.focalN))
        AN = args.focalN

    if not args.keepMono and (AC == 0 or AC == AN):
        continue

    focal_counts = [0, 0, 0, 0]
    focal_counts[NUCS.index(ref)] = AN - AC
    focal_counts[NUCS.index(alt)] = AC

    m = max(focal_counts)
    tied = [i for i, v in enumerate(focal_counts) if v == m]
    if len(tied) > 1:
        MAJOR, MINOR = NUCS[tied[0]], NUCS[tied[1]]
    else:
        MAJOR = NUCS[tied[0]]
        MINOR = alt if MAJOR == ref else ref

    row = [chrom, pos, ref, alt, MAJOR, MINOR,
           ",".join(map(str, focal_counts))]
    out.write("\t".join(row) + "\n")
    n_written += 1

inVCF.close()
out.close()
print(f"Wrote {n_written} sites to {args.outfile}")