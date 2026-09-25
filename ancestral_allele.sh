#!/bin/bash

function set_up_intel_env {

CONDA_SUBDIR=osx-64 conda create -y -n est-sfs python=3.9
#activate env
source $HOME/miniforge3/bin/activate est-sfs
echo $CONDA_DEFAULT_ENV #double check the env is correct

#permanently locks that particular Conda environment to only accept osx-64 (Intel) packages
conda config --env --set subdir osx-64
#install
conda install -c bioconda est-sfs repeatmasker -y
conda install -c bioconda numpy bcftools bwa bedtools -y
#MANUAL
}

function set_up_cactus_env {

conda create -n cactus-env -c conda-forge -c bioconda python=3.11 cactus -y
source $HOME/miniforge3/bin/activate cactus-env
pip install toil

}

#RUN CACTUS WGA: ALIGN OUTGROUPS TO REFERENCE
function run_cactus {

REF=$1        # reference fasta (rufipogon.fa)
OG1=$2        # outgroup1 assembly fasta (barthii)
OG2=$3        # outgroup2 assembly fasta (glamaepatula)
OUTDIR=$4     # output directory
NAME1=$5      # short name for OG1 (e.g. W1943) - no spaces/special chars
NAME2=$6      # short name for OG2 (e.g. W1715)
REFNAME=$7    # short name for reference (e.g. sativa)

source $HOME/miniforge3/bin/activate cactus-env 
echo "USING ENV" $CONDA_DEFAULT_ENV

mkdir -p $OUTDIR/work

# --- 1. build the seqFile (Newick tree + genome paths) ---
# Tree: reference + two outgroups. Adjust topology to true relationships.
SEQFILE=$OUTDIR/cactus_seqfile.txt
cat > $SEQFILE <<EOF
((${REFNAME}:0.003,${NAME1}:0.003):0.05,${NAME2}:0.08);
${REFNAME}	$(readlink -f $REF)
${NAME1}	$(readlink -f $OG1)
${NAME2}	$(readlink -f $OG2)
EOF

echo "=== seqFile ==="; cat $SEQFILE; echo

# --- 2. run cactus ---
JOBSTORE=$OUTDIR/jobstore
HAL=$OUTDIR/rice_wga.hal

# remove stale jobstore if re-running
rm -rf $JOBSTORE

cactus \
    $JOBSTORE \
    $SEQFILE \
    $HAL \
    --workDir $OUTDIR/work \
    --maxCores 8 \
    --binariesMode local \
    --branchScale 1.0

echo "=== HAL output ==="; ls -lh $HAL

# --- 3. quick verification of the HAL ---
halStats $HAL

}

#EXTRACT COORDINATES FROM HAL
function get_vcf_coordinates {

VCF_FILE=$1
OUTDIR=$2     # output directory

source $HOME/miniforge3/bin/activate est-sfs
echo "USING ENV" $CONDA_DEFAULT_ENV

echo "Step 1: Extracting SNP positions from VCF ..."
# Create a BED file of SNP positions (0-based start, 1-based end for hal2maf --targets)
# VCF is 1-based POS. BED is 0-based start.
# We want a 1bp window: start = POS-1, end = POS
bcftools query -f '%CHROM\t%POS\n' "$VCF_FILE"  \
    |awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2,$1":"$2}' \
    |grep -v 'ChrUN'> $OUTDIR/ref_snps.bed

# Also get the Ref allele to verify later (optional but good for sanity)
bcftools query -f '%CHROM\t%POS\t%REF\n' "$VCF_FILE" > $OUTDIR/snps.ref.tmp.tsv
}

function get_hal_coordinates {

HAL=$1
BED=$2
REFNAME=$3    # short name for reference (e.g. sativa)
OGNAME=$4      #short name for OG1 (e.g. W1943) - no spaces/special chars
OUTDIR=$5     # output directory

source $HOME/miniforge3/bin/activate cactus-env
echo "USING ENV" $CONDA_DEFAULT_ENV

echo "Step 2: Extracting SNP positions on OG from HAL on ${OGNAME}..."
halLiftover --bedType 4 $HAL $REFNAME $BED $OGNAME $OUTDIR/snps_${OGNAME}.bed
}
function get_OG_bases {

BED=$1
OG_FA=$2
OGNAME=$3
OUTDIR=$4

source $HOME/miniforge3/bin/activate cactus-env
echo "USING ENV" $CONDA_DEFAULT_ENV

echo "Step 3: Extracting base on OG from HAL on ${OGNAME}..."
bedtools getfasta -fi $OG_FA \
                  -bed $BED \
                  -tab -nameOnly > $OUTDIR/snps_${OGNAME}.tsv
}

#CONVERT TO EST-SFS INPUT (OG PART)
function sfs_OG {

REF_FILE=$1
OG_FILE=$2
OGNAME=$3
OUTDIR=$4

echo "Step 3: Make allele freq of ${OGNAME}..."
# Use awk to build a hash of ref positions and merge with outgroup bases
awk -v OFS='\t' '
    # Load outgroup file first (NR==FNR reads first file)
    NR==FNR { 
        og[$1] = $2 #   Store a hash: og[chrom:pos] = base
        next
    }
    # Process reference file
    {
        key = $1":"$2          # Build a lookup key from ref: "chrom:pos"e.g. "1:1118"
        base = og[key]         # look up outgroup base
        
        if (base == "A")      acgt = "1,0,0,0"
            else if (base == "C") acgt = "0,1,0,0"
            else if (base == "G") acgt = "0,0,1,0"
            else if (base == "T") acgt = "0,0,0,1"
            else { base = "NA"; acgt = "0,0,0,0" }

            ogcoord = (key in og_pos) ? og_pos[key] : "NA"

            # Columns:
            # 1 chrom_ref  2 pos_ref  3 base_ref
            # 4 og_coord   5 og_base
            # 6 focal-style ACGT vector (outgroup)
            print $1, $2, $3, ogcoord, base, acgt
    }
' $OG_FILE $REF_FILE > $OUTDIR/allele_count_$OGNAME.txt

echo "Done. Output written to $OUTDIR/allele_count_$OGNAME.txt"
echo "Total positions: $(wc -l < $OUTDIR/allele_count_$OGNAME.txt)"

}


#CONVERT PLINK TO VCF
function plink_to_vcf {

INDIR=$1
OUTDIR=$2
PREFIX1=$3
PREFIX2=$4

source $HOME/miniforge3/bin/activate est-sfs
echo "USING ENV" $CONDA_DEFAULT_ENV

for i in $(seq 1 2); do
    plink2 \
        --pfile $INDIR/$PREFIX1 \
        --chr $i \
        --export vcf bgz id-paste=iid \
        --out $OUTDIR/${PREFIX2}.chr${i}
    
    # strip INFO right after export (focal is already GT-only in FORMAT)
    bcftools annotate -x INFO --threads 4 \
        $OUTDIR/${PREFIX2}.chr${i}.vcf.gz \
        -Oz -o $OUTDIR/${PREFIX2}.chr${i}.noinfo.vcf.gz
    #RID OF THE FILE WITH INFO
    mv $OUTDIR/${PREFIX2}.chr${i}.noinfo.vcf.gz $OUTDIR/${PREFIX2}.chr${i}.vcf.gz
done
}

#MAKE ALLELE COUNT INPUT FOR est-sfs
function count_allele {

INFILE=$1
OUTDIR=$2
NAME=$3

source $HOME/miniforge3/bin/activate est-sfs
echo "USING ENV" $CONDA_DEFAULT_ENV
echo "CREATING est-sfs INPUT ALLELE COUNT FILES"

#make focal sample list
#bcftools query -l $INFILE |grep -v -E '^(W1943|W1715)$' > focal.txt 

python3 prep_estsfs.py \
    --infile $INFILE \
    --focalSampleList focal.txt \
    --outfile $OUTDIR/allele_count.${NAME}.tsv \
    --focalN 100 \
    --keepMono \
    --seed 42

cat $OUTDIR/allele_count.${NAME}.tsv|cut -f 7- > $OUTDIR/estsfs_data.${NAME}.txt
}


function run_est_sfs {

DATA=$1
CONFIG=$2
OUTDIR=$3
NAME=$4

source $HOME/miniforge3/bin/activate est-sfs
echo "USING ENV" $CONDA_DEFAULT_ENV

est-sfs $CONFIG $DATA seedfile.txt \
        ${OUTDIR}/${NAME}_sfs_out.txt \
        ${OUTDIR}/${NAME}_pvalue_out.txt

}

#MASKING
function run_repeatmasker {

FASTA=$1

RepeatMasker $FASTA \
    -species ''Oryza'' \
    -xsmall \
    -pa 8 \
    -gff
}


#---------------------execute------------------------

function main {  

OUTGRP=/Users/zchen13/data/outgroups
W1943=$OUTGRP/w1943-PCRfree_Round82_Lane2_sorted.uq_removeIndels.hardFiltered.vcf.gz
W1715=$OUTGRP/w1715_sorted.uq_removeIndels.hardFiltered.vcf.gz
SAMPLE=/Users/zchen13/data/oryza_19K_RGP
IRGSP=/Users/zchen13/data/reference/IRGSP-1.0_chr.renamed.fa
BARTHII=/Users/zchen13/data/reference/O_barthii.chr.fa
GLAMPTL=/Users/zchen13/data/reference/O_glamaepatula.chr.fa
WORKDIR=/Users/zchen13/SCRATCH/ancestral_allele
#
RESULT1=$WORKDIR/01_est-sfs
RESULT2=$WORKDIR/02_

#OUTPUT SHORTCUTS
VCF_ALL=$RESULT1/genome1.focal.chr1-12.vcf.gz

#est-sfs config files
JC=/Users/zchen13/scripts/ancestral-allele/config-JC.txt
KIMURA=/Users/zchen13/scripts/ancestral-alleleconfig-kimura.txt
R6=/Users/zchen13/scripts/ancestral-allele/config-rate6.txt

#set_up_intel_env

function setup_dir {

mkdir -p $WORKDIR
mkdir -p $RESULT1
#mkdir -p $RESULT2
#mkdir -p $RESULT3
#mkdir -p $RESULT5
#mkdir -p $RESULT5

}
setup_dir

#CONVERT PLINK TO VCF
#plink_to_vcf $SAMPLE $RESULT1 genome1.biallelic.base.genomewide.SNPs.withID.PlinkFormat genome1.focal

#MAKE SURE ONLY INCLUDE BIALLILIC SNPs, rename chr and extract chr 1-12
#biallelic $RESULT1/genome1.focal.vcf.gz $RESULT1/genome1.focal.biallelic.vcf.gz

#WGA WITH CACTUS FOR OUTLIER ALLELE EXTRACTION
#run_repeatmasker $GLAMPTL
#run_cactus $IRGSP $BARTHII $GLAMPTL $RESULT1 Obrth Oglmptl Ostv

#OUTGROUP ALLELE EXTRACTION
#get_vcf_coordinates $VCF_ALL $RESULT1
#get_hal_coordinates $RESULT1/rice_wga.hal $RESULT1/ref_snps.bed Ostv Obrth $RESULT1
#get_hal_coordinates $RESULT1/rice_wga.hal $RESULT1/ref_snps.bed Ostv Oglmptl $RESULT1  
#get_OG_bases $RESULT1/snps_Obrth.bed $BARTHII Obrth $RESULT1
#get_OG_bases $RESULT1/snps_Oglmptl.bed $GLAMPTL Oglmptl $RESULT1
#sfs_OG $RESULT1/snps_IRGSP.tsv $RESULT1/snps_Obrth.tsv Obrth $RESULT1
sfs_OG $RESULT1/snps_IRGSP.tsv $RESULT1/snps_Oglmptl.tsv Oglmptl $RESULT1

#GENERATE EST-SFS INPUT FILES
#count_allele $VCF_ALL $RESULT1 19k_samples

#RUN EST-SFS
#run_est_sfs $RESULT1/estsfs_data.chr1.tsv $R6 $RESULT1 chr1
#run_est_sfs chr1_first100.tsv $R6 $RESULT1 chr1
}


main