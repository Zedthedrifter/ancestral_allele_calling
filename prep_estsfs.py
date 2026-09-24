#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Prepare an input file for Keightley & Jackson (2018) est-sfs.

For every biallelic SNP in a (optionally gzipped) VCF, writes one row:
  chrom  pos  ref  alt  MAJOR  MINOR  focal[A,C,G,T]  og1[A,C,G,T]  og2[A,C,G,T]

The last three columns (7-9), stripped and space-separated, ARE the est-sfs
data file. The leading columns are kept for bookkeeping / mapping est-sfs
output back to genomic coordinates and to your major/minor alleles.

Usage:
  python3 prep_estsfs.py \
      --infile in.vcf.gz \
      --focalSampleList focal.txt \
      --outgroup1SampleList og1.txt \
      --outgroup2SampleList og2.txt \
      --outfile out.tsv \
      --focalN 100 --seed 42          # <- constant ingroup size for est-sfs
"""
import gzip
import argparse
import re
import numpy as np

# ----------------------------------------------------------------------
# Arguments
# ----------------------------------------------------------------------
p = argparse.ArgumentParser(description="Generate est-sfs input from a VCF")
p.add_argument("--infile", required=True)
p.add_argument("--focalSampleList", required=True)
p.add_argument("--outgroup1SampleList", required=True)
p.add_argument("--outgroup2SampleList", required=True)
p.add_argument("--outfile", required=True)
# est-sfs needs a CONSTANT ingroup sample size. If set, every site is randomly
# down-sampled (without replacement) to this many alleles; sites with fewer
# called alleles are dropped. Strongly recommended for real data.
p.add_argument("--focalN", type=int, default=None,
               help="target number of focal ALLELES (constant n for est-sfs)")
p.add_argument("--seed", type=int, default=42,
               help="RNG seed for down-sampling and 0.5-freq outgroup coin flips")
p.add_argument("--keepMono", action="store_true",
               help="also write sites that are monomorphic in the focal group")
args = p.parse_args()

rng = np.random.default_rng(args.seed)   # single seeded RNG -> reproducible

# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def smart_open(path):
    """Open plain or gzipped text transparently (fixes B5)."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, "rt")

def read_list(path):
    with open(path) as fh:
        return fh.read().splitlines()

GT_SPLIT = re.compile(r"[/|]")   # split on both phased '|' and unphased '/'

def count_ref_alt(calls):
    """
    Count called REF(0) and ALT(1) alleles across a list of GT fields.
    Robust to phasing, half-calls, haploid calls, and any missing form,
    because it inspects individual alleles rather than matching whole-GT
    strings (fixes B2, B3). Alleles that are '.' or index >=2 are ignored;
    the latter should not occur once the VCF is split to biallelic SNPs.
    Returns (AN = called alleles, AC = alt allele count).
    """
    an = ac = 0
    for c in calls:
        gt = c.partition(":")[0]          # genotype is the first ':'-field
        for a in GT_SPLIT.split(gt):
            if a == "0":
                an += 1
            elif a == "1":
                an += 1
                ac += 1
            # '.', '2', etc. -> not counted
    return an, ac

NUCS = ["A", "C", "G", "T"]               # order REQUIRED by est-sfs

# ----------------------------------------------------------------------
# Load sample lists and check they are mutually exclusive
# ----------------------------------------------------------------------
focal = read_list(args.focalSampleList)
og1   = read_list(args.outgroup1SampleList)
og2   = read_list(args.outgroup2SampleList)

for a, b, name in [(focal, og1, "focal/og1"),
                   (focal, og2, "focal/og2"),
                   (og1,   og2, "og1/og2")]:
    if set(a) & set(b):
        raise SystemExit(f"Sample lists {name} overlap!")

# ----------------------------------------------------------------------
# Read the VCF header once, build name -> column-index map (fixes B7)
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
missing = [s for grp in (focal, og1, og2) for s in grp if s not in sample_idx]
if missing:
    raise SystemExit(f"Samples not in VCF: {missing[:10]} ...")

# Column indices (into the genotype fields, i.e. line[9:]) computed ONCE
focal_cols = [sample_idx[s] for s in focal]
og1_cols   = [sample_idx[s] for s in og1]
og2_cols   = [sample_idx[s] for s in og2]

# ----------------------------------------------------------------------
# Main loop
# ----------------------------------------------------------------------
out = open(args.outfile, "w")            # plain text; bgzip later if wanted
n_written = 0

for line in inVCF:
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    chrom, pos, ref, alt = f[0], f[1], f[3], f[4]

    # Keep only biallelic SNPs with single-nucleotide REF/ALT (handles B6)
    if ref not in NUCS or alt not in NUCS:
        continue

    geno = f[9:]

    # ---- focal counts ----
    focal_calls = [geno[i] for i in focal_cols]
    AN, AC = count_ref_alt(focal_calls)
    if AN == 0:
        continue

    # ---- enforce constant ingroup size for est-sfs (the key fix) ----
    if args.focalN is not None:
        if AN < args.focalN:
            continue                                  # too little data -> drop
        # hypergeometric draw = alt count when sampling focalN of AN alleles
        AC = int(rng.hypergeometric(AC, AN - AC, args.focalN))
        AN = args.focalN

    #if not args.keepMono and (AC == 0 or AC == AN):
    #    continue                                      # skip focal-monomorphic

    focal_counts = [0, 0, 0, 0]
    focal_counts[NUCS.index(ref)] = AN - AC
    focal_counts[NUCS.index(alt)] = AC

    # ---- major / minor using est-sfs's convention (ties -> alphabetical) ----
    m = max(focal_counts)
    tied = [i for i, v in enumerate(focal_counts) if v == m]
    if len(tied) > 1:                                 # equal freq
        MAJOR, MINOR = NUCS[tied[0]], NUCS[tied[1]]
    else:
        MAJOR = NUCS[tied[0]]
        MINOR = alt if MAJOR == ref else ref

    # ---- outgroups: place a single '1' on the major allele ----
    def og_vector(cols):
        calls = [geno[i] for i in cols]
        an, ac = count_ref_alt(calls)
        v = [0, 0, 0, 0]
        if an == 0:
            return v                                   # no info -> all zeros
        frac = ac / an
        if frac == 0.5:
            major = ref if rng.random() < 0.5 else alt # seeded coin flip
        else:
            major = alt if frac > 0.5 else ref
        v[NUCS.index(major)] = 1
        return v

    og1_v = og_vector(og1_cols)
    og2_v = og_vector(og2_cols)

    # ---- write ----
    row = [chrom, pos, ref, alt, MAJOR, MINOR,
           ",".join(map(str, focal_counts)),
           ",".join(map(str, og1_v)),
           ",".join(map(str, og2_v))]
    out.write("\t".join(row) + "\n")
    n_written += 1

inVCF.close()
out.close()
print(f"Wrote {n_written} sites to {args.outfile}")

