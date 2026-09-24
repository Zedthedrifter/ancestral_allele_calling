# -*- coding: utf-8 -*-
"""
Spyder Editor

This script will prepare an infile for Keightley's ext SFS method for determining the ancestral allele

requires 2 outgroups (for now)

and should output chrom and position for now (use awk to get rid of prior to ext sfs?)

missing data is okay

Note: I am still writing out ambiguous file because it's a handy output for parsing those sites downstream.
script usage:
python3 script --infile [] \
    --focalSampleList [] \
          --outgroup1SampleList [] \
               --outgroup2SampleList [] \
                   --outfile []
"""
import sys
import gzip
import argparse
import random # for coin flips

##### set args and sample lists
from typing import Any

parser = argparse.ArgumentParser(description='Script to generate input for extsfs from a vcf file')
parser.add_argument("--infile",required=True,help="path to vcf file containing all your sequences.")
parser.add_argument("--focalSampleList",required=True,help="path to file containing IDs of samples that are the focal species. Should be unrelated and from the same population.")
parser.add_argument("--outgroup1SampleList",required=True,help="path to file containing IDs of samples that are the closer species to your focal species of the two outgroups")
parser.add_argument("--outgroup2SampleList",required=True,help="path to file containing IDs of samples that are the furthest species to your focal species of the two outgroups")
parser.add_argument("--outfile",required=True,help="path to output gzipped vcf file")

args = parser.parse_args()

infile=str(args.infile)

focalSamplesFile= open(str(args.focalSampleList), "r")
focalSamples = focalSamplesFile.read().splitlines()
focalSamplesFile.close()

outgroup1SamplesFile= open(str(args.outgroup1SampleList), "r")
outgroup1Samples = outgroup1SamplesFile.read().splitlines()
outgroup1SamplesFile.close()

outgroup2SamplesFile= open(str(args.outgroup2SampleList), "r")
outgroup2Samples = outgroup2SamplesFile.read().splitlines()
outgroup1SamplesFile.close()

outfilename=str(args.outfile)

#### make sure sample sets are mutually exclusive
if sum( [x in focalSamples for x in outgroup1Samples])>0:
    print("your focal and outgroup1 lists are not exclusive!",flush=True)
    exit(1)
if sum( [x in focalSamples for x in outgroup2Samples])>0:
    print("your focal and outgroup2 lists are not exclusive!",flush=True)
    exit(1)
if sum( [x in outgroup1Samples for x in outgroup2Samples])>0:
    print("your outgroup1 and outgroup2 lists are not exclusive!",flush=True)
    exit(1)


inVCF=gzip.open(infile, 'rt')
outfile = open(outfilename, 'w') # wanted to write out to gzip file but that doesn't work with tabix. need bgzipped file annoyingly.
# don't have a header -- messes stuff up downstream
# get samples of vcf
samples=[]
for line in inVCF:
    if line.startswith('##'):
        pass
    else:
        for i in line.split()[9:]: samples.append(i)
        break
inVCF.seek(0)


# make sure all your samples are in the list:
if sum([ x in samples for x in focalSamples])==len(focalSamples):
    print("all focal samples are in the vcf",flush=True)
else:
    print("some reference samples are missing from vcf!",flush=True)
    exit(1)
if sum([ x in samples for x in outgroup1Samples])==len(outgroup1Samples):
    print("all outgroup1 samples are in the vcf",flush=True)
else:
    print("some outgroup1 samples are missing from vcf!",flush=True)
    exit(1)
if sum([ x in samples for x in outgroup2Samples])==len(outgroup2Samples):
    print("all outgroup2 samples are in the vcf",flush=True)
else:
    print("some outgroup2 samples are missing from vcf!",flush=True)
    exit(1)

# make sure there are no extraneous individuals:
if sum([ x not in focalSamples+outgroup1Samples+outgroup2Samples for x in samples])>0:
    print("There are individuals present in your vcf that aren't in your reference/sister/outgroup lists -- that should be ok but just letting you know!")
    #exit(1)
else:
    print("there are no extraneous samples in your vcf")

print(str("focal sample:"+str(focalSamples)),flush=True)
print(str("outgroup1 samples:"+str(outgroup1Samples)),flush=True)
print(str("outgroup2 samples:"+str(outgroup2Samples)),flush=True)


# Define nucleotides
nucs=['A','C','G','T'] # this is the order they must be in for ext sfs

######## go through vcf and estimate ancestral alleles ##########
for line0 in inVCF:
    if line0.startswith('#'):
        continue
    AA=None # reset each time
    line=line0.strip().split('\t') # this splits line by tabs
    mychr=line[0]
    mypos=line[1]
    myref=line[3]
    myalt=line[4]
    if myalt not in nucs:
        #print(str(mychr)+":",str(mypos)+" alt allele: "+str(myalt)+" is not in "+str(nucs)+", skipping it!")
        # skip to next site
        continue
    if myref not in nucs:
        #print("ref allele: " + str(myref) + " is not in " + str(nucs) + ", skipping it!")
        # skip to next site
        continue
    mygenoinfo=line[9:]
    # Initialize counts in order of nucs
    focal_counts = [0, 0, 0, 0]
    outgroup1_counts = [0, 0, 0, 0]
    outgroup2_counts = [0, 0, 0, 0]
    # get genotype info
    allCalls=[i.split(":")[0] for i in mygenoinfo] # get genotype calls
    callDict = dict(zip(samples,allCalls)) # associate with sample names
    # pull out samples in each category
    focalSpeciesCalls=[callDict[x] for x in focalSamples]
    outgroup1SpeciesCalls=[callDict[x] for x in outgroup1Samples]
    outgroup2SpeciesCalls=[callDict[x] for x in outgroup2Samples]
    # need AN and AC for focal species
    focalAN=2*(len(focalSpeciesCalls) - focalSpeciesCalls.count("./.")) # 2*(totalNonMissingGTcalls)
    focalHet=focalSpeciesCalls.count("0/1")+focalSpeciesCalls.count("1/0")+focalSpeciesCalls.count("0|1")+focalSpeciesCalls.count("1|0")
    focalHomAlt=focalSpeciesCalls.count("1/1")+focalSpeciesCalls.count("1|1")
    focalAC=focalHet+2*focalHomAlt
    # if there are no non-missing calls skip to the next site in the vcf
    if focalAN == 0:
        continue
    # if focalAN == focalAC then the site is monomorphic for the focal species
    # index based on the nucslist -- get the index of the ref allele
    # soo ref allele counts AN- AC and alt allele count is AC (works even if AC = 0)
    #ref allele:
    focal_counts[nucs.index(myref)] = int(focalAN) - int(focalAC)
    # alt allele:
    focal_counts[nucs.index(myalt)] = int(focalAC)
    # check:
    if sum(focal_counts) != focalAN:
        print("something is weird -- nuc count sum doesn't equal focal AN",flush=True)
        exit(1)

    # for downstream use, want to get the major and minor allele in the focal species. Note that if they have the exact same frequency, ext-sfs treats the major allele as the one that comes first alphabetically (A,C,G,T)
    # Get the position of the most frequent allele in the focal group -- give the position 0,1,2,3
    max_pos=[i for i,j in enumerate(focal_counts) if j==max(focal_counts)] # this works because enumerate gives paired values of location in interation through list,item so you are wanting the first entry (location) which is i, only if j (the second entry) equals the max value. You have to do it this way because otherwise you might miss sites where the major/minor alleles have equal frequencies
    # if there is more than 1 max position, means major freq = minor freq. so extsfs chooses the alphabetically first one as the major allele, so I will do so here as well: the positions are in alphabetical order in max_pos so the first one is major, second is minor
    if len(max_pos) > 1:
        MAJOR = nucs[max_pos[0]]
        MINOR = nucs[max_pos[1]]
    else:
        MAJOR=nucs[max_pos[0]]
        if MAJOR==myref:
            MINOR=myalt
        elif MAJOR==myalt:
            MINOR=myref
        else:
            print("major allele not == ref or alt, what's going on?",flush=True)
            exit(1)

    ### for outgroups, need to get major allele:
    # outgroup 1
    outgroup1AN=2*(len(outgroup1SpeciesCalls) - outgroup1SpeciesCalls.count("./.")) # 2*(totalNonMissingGTcalls)
    outgroup1Het=outgroup1SpeciesCalls.count("0/1")+outgroup1SpeciesCalls.count("1/0")+outgroup1SpeciesCalls.count("0|1")+outgroup1SpeciesCalls.count("1|0")
    outgroup1HomAlt=outgroup1SpeciesCalls.count("1/1")+outgroup1SpeciesCalls.count("1|1")
    outgroup1AC=outgroup1Het+2*outgroup1HomAlt
    # figure out which one is the major allele
    # if there are no calls for the outgroup you put 0,0,0,0
    if outgroup1AN == 0:
        outgroup1_counts=[0,0,0,0]
    # otherwise get
    elif float(outgroup1AC) / float(outgroup1AN) == 0.5:
        # equal amounts so flip a coin
        outgroup1majorAllele=random.choice([myref, myalt]) # choose randomly between ref and alt to be major allele
        outgroup1_counts[nucs.index(outgroup1majorAllele)] = 1
    elif float(outgroup1AC) / float(outgroup1AN) < 0.5:
        # alt allele is minor allele so ref is major allele
        outgroup1majorAllele = myref
        outgroup1_counts[nucs.index(outgroup1majorAllele)] = 1
    elif float(outgroup1AC) / float(outgroup1AN) > 0.5:
        # alternate allele is major allele and ref is minor allele
        outgroup1majorAllele = myalt
        outgroup1_counts[nucs.index(outgroup1majorAllele)] = 1
    else:
        print("something is strange with your outgroup1 at this site",flush=True)
        exit(1)

    # outgroup 2:

    outgroup2AN=2*(len(outgroup2SpeciesCalls) - outgroup2SpeciesCalls.count("./.")) # 2*(totalNonMissingGTcalls)
    outgroup2Het=outgroup2SpeciesCalls.count("0/1")+outgroup2SpeciesCalls.count("1/0")+outgroup2SpeciesCalls.count("0|1")+outgroup2SpeciesCalls.count("1|0")
    outgroup2HomAlt=outgroup2SpeciesCalls.count("1/1")+outgroup2SpeciesCalls.count("1|1")
    outgroup2AC=outgroup2Het+2*outgroup2HomAlt

    # figure out which one is the major allele
    # if there are no calls for the outgroup you put 0,0,0,0
    if outgroup2AN == 0:
        outgroup2_counts=[0,0,0,0]
    # otherwise get
    elif float(outgroup2AC) / float(outgroup2AN) == 0.5:
        # equal amounts so flip a coin
        outgroup2majorAllele=random.choice([myref, myalt]) # choose randomly between ref and alt to be major allele
        outgroup2_counts[nucs.index(outgroup2majorAllele)] = 1
    elif float(outgroup2AC) / float(outgroup2AN) < 0.5:
        # alt allele is minor allele so ref is major allele
        outgroup2majorAllele = myref
        outgroup2_counts[nucs.index(outgroup2majorAllele)] = 1
    elif float(outgroup2AC) / float(outgroup2AN) > 0.5:
        # alternate allele is major allele and ref is minor allele
        outgroup2majorAllele = myalt
        outgroup2_counts[nucs.index(outgroup2majorAllele)] = 1
    else:
        print("something is strange with your outgroup2 at this site",flush=True)
        exit(1)

    # okay time to write out: (adding major and minor)
    output = [str(mychr), str(mypos), str(myref),str(myalt),str(MAJOR),str(MINOR),",".join(str(count) for count in focal_counts), ",".join(str(count) for count in outgroup1_counts), ",".join(str(count) for count in outgroup2_counts)]
    #print("\t".join(output))
    outfile.write("\t".join(output)+"\n")

inVCF.close()
outfile.close()
