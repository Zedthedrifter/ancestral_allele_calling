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
NAME1=$5      # short name for OG1 (e.g. Obrthii) - no spaces/special chars
NAME2=$6      # short name for OG2 (e.g. Oglmptl)
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

function filter_est-sfs_inputs {

#FILTER: MAKE SURE THE SITES ARE IN SAME ORDER FOR EST-SFS
awk -F'\t' 'NR==FNR { keep[$1":"$2]=1; next } ($1":"$2) in keep' \
allele_count.19k_samples.tsv allele_count_Obrth.txt > allele_count_filtered_Obrth.txt
awk -F'\t' 'NR==FNR { keep[$1":"$2]=1; next } ($1":"$2) in keep' \
allele_count.19k_samples.tsv allele_count_Oglmptl.txt > allele_count_filtered_Oglmptl.txt

#COMPILE: IN EST-SFS INPUT FORMAT

}

function compile_est-sfs_inputs {

FOCAL=$1
OG1=$2
OG2=$3
OUTDIR=$4

echo "USING ENV" $CONDA_DEFAULT_ENV
echo "COMPILING INPUT ALLELE COUNT FILES"

cat $FOCAL|grep -v 'ChrUN'|cut -f 7 > $OUTDIR/tmp1.tsv
cut -f 6 $OG1 > $OUTDIR/tmp2.tsv
cut -f 6 $OG2 > $OUTDIR/tmp3.tsv
wc -l $OUTDIR/tmp*
paste $OUTDIR/tmp1.tsv $OUTDIR/tmp2.tsv $OUTDIR/tmp3.tsv > $OUTDIR/est_sfs_input.tsv
rm -f $OUTDIR/tmp*
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

#RUN IN BATCH TO AVOID MEMORY ISSUE
function run_est_sfs_chunked {

DATA=$1        # full est-sfs input file
CONFIG=$2      # config file (e.g. config_r6.txt)
OUTDIR=$3      # output directory
NAME=$4        # prefix for output names (e.g. rice)
CHUNK=$5       # lines per chunk (default: 2000000)

# Default chunk size if not provided
CHUNK=${CHUNK:-2000000}

source $HOME/miniforge3/bin/activate est-sfs
echo "USING ENV" $CONDA_DEFAULT_ENV

# Create a working directory for chunks
CHUNKDIR=$OUTDIR/chunks_${NAME}
mkdir -p $CHUNKDIR
LOG=$CHUNKDIR/errors.log

echo "Splitting $DATA into chunks of $CHUNK lines..."
split -l $CHUNK "$DATA" "${CHUNKDIR}/chunk_"

# Loop over each chunk
for chunk in ${CHUNKDIR}/chunk_*; do
    # Extract chunk number (e.g. chunk_00.txt -> 00)
    chunk_id=$(basename "$chunk" | sed 's/chunk_//')
    
    echo "Running est-sfs on chunk ${chunk_id}..."
    
    est-sfs $CONFIG $chunk seedfile.txt \
        ${CHUNKDIR}/${NAME}_sfs_out_${chunk_id}.txt \
        ${CHUNKDIR}/${NAME}_pvalue_out_${chunk_id}.txt
    
    if [ $? -ne 0 ]; then
        echo "ERROR: est-sfs failed on chunk ${chunk_id}"
        echo "ERROR: est-sfs failed on chunk ${chunk_id}" >> $LOG
    fi
done

echo "All chunks completed. Outputs in ${CHUNKDIR}/"

# Optional: concatenate all outputs into one file
echo "Concatenating outputs..."
for f in ${CHUNKDIR}/${NAME}_sfs_out_*.txt; do
    tr -d '\n' < "$f"
    printf ','
done > ${OUTDIR}/${NAME}_sfs_out_all.txt
#manually remove the last comma...

# PVALUE: strip first 8 header lines from each chunk, then concatenate
# (only keep the header from the FIRST chunk)
first=1
for f in ${CHUNKDIR}/${NAME}_pvalue_out_*.txt; do
    if [ $first -eq 1 ]; then
        cat "$f" >> ${OUTDIR}/${NAME}_pvalue_out_all.txt
        first=0
    else
        tail -n +9 "$f" >> ${OUTDIR}/${NAME}_pvalue_out_all.txt
    fi
done

echo "Done. Combined outputs:"
echo "  ${OUTDIR}/${NAME}_sfs_out_all.txt"
echo "  ${OUTDIR}/${NAME}_pvalue_out_all.txt"

# Optional: clean up chunks (uncomment to enable)
# rm -rf ${CHUNKDIR}

}

function extract_ancestral_allele {

THRESHOLD1=$1    # certainty for max_p (base identity)
THRESHOLD2=$2    # certainty for P-major-ancestral
PVALUE=$3        # est-sfs pvalue output
ALLELECOUNT=$4   # allele count file (focal, og1, og2)
OUTDIR=$5
OUTFILE=$6

echo "USING ENV" $CONDA_DEFAULT_ENV "to extract ancestral allele"
echo "THRESHOLD1 (max_p) =" $THRESHOLD1
echo "THRESHOLD2 (major_p) =" $THRESHOLD2

# ---- 1. Strip header lines from the pvalue file ----
PVALUE_DATA=${OUTDIR}/pvalue_data.tmp
awk '!/^0 /' $PVALUE > $PVALUE_DATA

# ---- 2. Confirm line counts match ----
N_PVALUE=$(wc -l < $PVALUE_DATA)
N_ALLELECOUNT=$(wc -l < $ALLELECOUNT)

if [ "$N_PVALUE" -ne "$N_ALLELECOUNT" ]; then
    echo "ERROR: line count mismatch"
    echo "  pvalue data lines:    $N_PVALUE"
    echo "  allele count lines:   $N_ALLELECOUNT"
    rm -f $PVALUE_DATA
    return 1
fi

echo "Line counts match: $N_PVALUE sites"

# ---- 3. Paste the two files side by side ----
JOINED=${OUTDIR}/joined.tmp
paste $ALLELECOUNT $PVALUE_DATA > $JOINED
echo "computing ancestral allele"
# ---- 4. Parse the joined file in a single pass ----
awk -v THRESHOLD1="$THRESHOLD1" -v THRESHOLD2="$THRESHOLD2" '
    {
        # --- Fields from allele count file (columns 1-3) ---
        chr          = $1
        site         = $2
        focal_counts = $3
        og1          = $4
        og2          = $5

        # --- Fields from pvalue file (columns 6-12) ---
        xxxx  = $6
        code  = $7
        p_maj = $8
        p_A   = $9
        p_C   = $10
        p_G   = $11
        p_T   = $12

        # --- Focal counts from the comma-delimited string ---
        split(focal_counts, f, ",")
        fa = f[1]; fc = f[2]; fg = f[3]; ft = f[4]

        # --- Most likely ancestral base from pvalue ---
        max_p = p_A; anc = "A"
        if (p_C > max_p) { max_p = p_C; anc = "C" }
        if (p_G > max_p) { max_p = p_G; anc = "G" }
        if (p_T > max_p) { max_p = p_T; anc = "T" }

        # --- Major allele from focal counts ---
        maj_p = fa; maj = "A"
        if (fc > maj_p) { maj_p = fc; maj = "C" }
        if (fg > maj_p) { maj_p = fg; maj = "G" }
        if (ft > maj_p) { maj_p = ft; maj = "T" }

        # --- Decision rule ---
        if (max_p > THRESHOLD1) {
            aa = anc
        } else if (p_maj > THRESHOLD2) {
            aa = maj
        } else {
            aa = "NA"
        }

        # --- Output ---
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", \
            chr, site, aa, maj, p_maj, \
            p_A, p_C, p_G, p_T, \
            focal_counts, og1, og2
    }
' $JOINED > ${OUTDIR}/${OUTFILE}

# ---- 5. Clean up temporaries ----
rm -f $PVALUE_DATA $JOINED

echo "Output written to ${OUTDIR}/${OUTFILE}"
}

#PLOT uSFS

function plot_uSFS {

AA_FILE=$1
NAME=$2
OUTDIR=$3

source $HOME/miniforge3/bin/activate est-sfs
echo "USING ENV" $CONDA_DEFAULT_ENV

python3 plot_uSFS.py $AA_FILE $NAME $OUTDIR
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
MRDONLS=/Users/zchen13/data/reference/O_meridionalis.chr.fa
WORKDIR=/Users/zchen13/SCRATCH/ancestral_allele
#
RESULT1=$WORKDIR/01_est-sfs
RESULT2=$WORKDIR/02_est-sfs

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
mkdir -p $RESULT2
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
#run_cactus $IRGSP $BARTHII $MRDONLS $RESULT2 Obrth Omrdonls Ostv

#OUTGROUP ALLELE EXTRACTION
#get_vcf_coordinates $VCF_ALL $RESULT1
#get_hal_coordinates $RESULT1/rice_wga.hal $RESULT1/ref_snps.bed Ostv Obrth $RESULT1
#get_hal_coordinates $RESULT1/rice_wga.hal $RESULT1/ref_snps.bed Ostv Oglmptl $RESULT1
#get_hal_coordinates $RESULT2/rice_wga.hal $RESULT1/ref_snps.bed Ostv Omrdonls $RESULT2   
#get_OG_bases $RESULT1/snps_Obrth.bed $BARTHII Obrth $RESULT1
#get_OG_bases $RESULT1/snps_Oglmptl.bed $GLAMPTL Oglmptl $RESULT1
#get_OG_bases $RESULT2/snps_Omrdonls.bed $MRDONLS Omrdonls $RESULT2
#sfs_OG $RESULT1/snps_IRGSP.tsv $RESULT1/snps_Obrth.tsv Obrth $RESULT1
#sfs_OG $RESULT1/snps_IRGSP.tsv $RESULT1/snps_Oglmptl.tsv Oglmptl $RESULT1
#sfs_OG $RESULT1/snps_IRGSP.tsv $RESULT2/snps_Omrdonls.tsv Omrdonls $RESULT2

#GENERATE EST-SFS INPUT FILE: SAMPLE
#count_allele $VCF_ALL $RESULT1 19k_samples
#TRIM THE OG FILES (CAUSE SOME OF THE SITES WERE FILTERED OUT IN THE 19K FILE DUE TO EXCESSIVE MISSING DATA)

#COMPILE
#compile_est-sfs_inputs $RESULT1/allele_count.19k_samples.tsv $RESULT1/allele_count_filtered_Obrth.txt $RESULT1/allele_count_filtered_Oglmptl.txt $RESULT1
#FILTER THE INPUTS TO MAKE SURE THE NUMBER OF SITES ARE THE SAME
#need to add this command to function later
#awk -F'\t' 'NR==FNR { keep[$1":"$2]=1; next } ($1":"$2) in keep' \
#    $RESULT1/allele_count.19k_samples.tsv $RESULT2/allele_count_Omrdonls.txt > $RESULT2/allele_count_filtered_Omrdonls.txt
#compile_est-sfs_inputs $RESULT1/allele_count.19k_samples.tsv $RESULT1/allele_count_filtered_Obrth.txt $RESULT2/allele_count_filtered_Omrdonls.txt $RESULT2

#RUN EST-SFS
#run_est_sfs $RESULT1/estsfs_data.chr1.tsv $R6 $RESULT1 chr1
#run_est_sfs $RESULT1/chr1_first100.tsv $R6 $RESULT1 chr1_test
#run_est_sfs_chunked $RESULT1/est_sfs_input.tsv $R6 $RESULT1 rice_chr1-12_est-sfs 200000
#run_est_sfs_chunked $RESULT2/est_sfs_input.tsv $R6 $RESULT2 rice_chr1-12_est-sfs 200000

#EXTRACT AA
#extract_ancestral_allele 0.95 0.85 $RESULT1/rice_chr1-12_est-sfs_pvalue_out_all.txt $RESULT1/allele_count_filtered_compiled.tsv $RESULT1 AA_Obrth-Oglmptl.tsv
#extract_ancestral_allele 0.95 0.85 $RESULT2/rice_chr1-12_est-sfs_pvalue_out_all.txt $RESULT2/allele_count_filtered_compiled.tsv $RESULT2 AA_Obrth-mrdonls.tsv

#PLOT usfs
#plot_uSFS $RESULT1/AA_Obrth-Oglmptl.tsv og_set1 $RESULT1
plot_uSFS $RESULT2/AA_Obrth-mrdonls.tsv  og_set2 $RESULT2

#PLOT DISAGREED AA
./disagreement_genome_wide_density_plot.py $RESULT2/disagreed_aa.tsv disagreement_density_plot

}


main