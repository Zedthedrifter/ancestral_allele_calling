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

#EXTRACT OUTGROUP ALLELE AT SNP SITES

function extract_wga_bases {

HAL_FILE=$1          # Path to rice_wga.hal
REF_GENOME=$2        # Name of reference genome in HAL (e.g., "sativa")
OG1_GENOME=$3        # Name of outgroup 1 in HAL (e.g., "barthii")
OG2_GENOME=$4        # Name of outgroup 2 in HAL (e.g., "glumaepatula")
VCF_FILE=$5          # Input VCF (focal SNPs)
OUT_FILE=$6          # Output TSV file

source $HOME/miniforge3/bin/activate cactus-env  # Activate env with hal tools
echo "USING ENV: $CONDA_DEFAULT_ENV"

# Check inputs
if [ ! -f "$HAL_FILE" ]; then echo "Error: HAL file not found"; return 1; fi
if [ ! -f "$VCF_FILE" ]; then echo "Error: VCF file not found"; return 1; fi

echo "Step 2: Extracting alignment blocks using hal2maf..."
# hal2maf extracts the alignment for specific targets.
# --refGenome: The genome we are querying coordinates in.
# --targetGenomes: The genomes we want bases from.
# --targets: The BED file of positions.
# --noAncestors: We only care about extant species.
# Output is MAF format (standard alignment format).

hal2maf "$HAL_FILE" /dev/stdout \
    --refGenome "$REF_GENOME" \
    --targetGenomes "$OG1_GENOME,$OG2_GENOME" \
    --targets snps.tmp.bed \
    --noAncestors \
    2>/dev/null | \
python3 -c "
import sys

# Configuration
ref_name = '$REF_GENOME'
og1_name = '$OG1_GENOME'
og2_name = '$OG2_GENOME'

current_chrom = None
current_pos = None
bases = {ref_name: None, og1_name: None, og2_name: None}

def flush_site():
    if current_chrom is not None:
        # Get bases, default to 'N' if missing (gap or no alignment)
        b_ref = bases[ref_name] if bases[ref_name] else 'N'
        b_og1 = bases[og1_name] if bases[og1_name] else 'N'
        b_og2 = bases[og2_name] if bases[og2_name] else 'N'
        print(f'{current_chrom}\t{current_pos}\t{b_ref}\t{b_og1}\t{b_og2}')

for line in sys.stdin:
    line = line.strip()
    if not line or line.startswith('#'):
        continue
    
    if line.startswith('a'):
        # New alignment block. If we were processing a site, flush it first
        # (Though for 1bp targets, usually one block per site or split)
        # Actually, MAF 'a' lines start a block. We process 's' lines within.
        # Reset for new block if needed, but let's handle logic in 's' lines.
        pass
        
    elif line.startswith('s'):
        parts = line.split()
        # Format: s <genome>.<chr> <start> <size> <strand> <srcSize> <seq>
        genome_part = parts[1]
        genome = genome_part.split('.')[0] # Handle genome.chrom naming
        chrom = parts[1].split('.')[1] if '.' in parts[1] else parts[1]
        start = int(parts[2])
        seq = parts[6].upper()
        
        # We are looking for 1bp alignments usually, but MAF might have context
        # Since we queried 1bp targets, the seq should be 1 char (or gap)
        # If the target was 1bp, the alignment row for that genome should correspond.
        
        # Map coordinates: MAF start is 0-based. VCF POS is 1-based.
        # If the block covers our target, extract the specific base.
        # For 1bp query, usually start == target_start.
        
        target_pos_1based = start + 1
        
        # Simple logic: if this row corresponds to the chromosome we are tracking
        # Note: hal2maf output might group multiple targets if close. 
        # But since we want specific POS, we rely on the order matching the BED or explicit coords.
        # Robust way: The BED order is preserved in hal2maf output usually.
        
        # Let's refine: We need to match the specific query position.
        # Since we passed a BED, the output blocks correspond to those regions.
        # We assume for 1bp targets, the 's' line seq is the base (or gap).
        
        base = seq[0] if len(seq) > 0 else 'N'
        if base == '-': base = 'N' # Treat gaps as missing
        
        # Store base for this genome for the current logical site
        # We need to group by position. 
        # Since hal2maf outputs blocks, and we queried 1bp, 
        # we can assume the first 's' lines of a block correspond to the target.
        
        # Better approach for script: 
        # Just collect bases for the current block. 
        # A block corresponds to one target region in our 1bp BED.
        if genome == ref_name:
            bases[ref_name] = base
            # Update current pos/chrom from reference row
            current_chrom = chrom
            current_pos = target_pos_1based
        elif genome == og1_name:
            bases[og1_name] = base
        elif genome == og2_name:
            bases[og2_name] = base
            
    elif line.startswith('//') or line == '':
        # End of block (sometimes marked) or empty line
        # Flush if we have data
        if current_chrom is not None and bases[ref_name]:
             flush_site()
             current_chrom = None
             bases = {ref_name: None, og1_name: None, og2_name: None}

# Final flush
if current_chrom is not None:
    flush_site()

" > "$OUT_FILE"

# Cleanup
#rm -f snps.tmp.bed snps.ref.tmp.tsv

echo "Done. Output saved to $OUT_FILE"
echo "First 10 lines:"
head "$OUT_FILE"

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

#RENAME OUTGROUPS 

function rename_outgrp {

echo W1715 > newname.txt
  bcftools view -m2 -M2 -v snps $W1715 -Ou \
| bcftools reheader -s newname.txt -o $OUTGRP/W1715.bcf
bcftools index $OUTGRP/W1715.bcf
#CHECK OUTPUT
echo "== W1715 =="
echo "orig samples:"; bcftools query -l $W1715 2>/dev/null | wc -l
echo "new samples :"; bcftools query -l $OUTGRP/W1715.bcf  2>/dev/null | wc -l
echo "orig total records:"; bcftools view -H $W1715 2>/dev/null | wc -l
echo "orig biallelic SNPs:"; bcftools view -H -m2 -M2 -v snps $W1715 2>/dev/null | wc -l
echo "new BCF records:"; bcftools view -H $OUTGRP/W1715.bcf 2>/dev/null | wc -l
bcftools query -l $OUTGRP/W1715.bcf 2>/dev/null

#now do W1943
echo W1943 > newname.txt
  bcftools view -m2 -M2 -v snps $W1943 -Ou \
| bcftools reheader -s newname.txt -o $OUTGRP/W1943.bcf
bcftools index $OUTGRP/W1943.bcf
#CHECK OUTPUT
echo "== W1715 =="
echo "orig samples:"; bcftools query -l $W1943 2>/dev/null | wc -l
echo "new samples :"; bcftools query -l $OUTGRP/W1943.bcf  2>/dev/null | wc -l
echo "orig total records:"; bcftools view -H $W1943 2>/dev/null | wc -l
echo "orig biallelic SNPs:"; bcftools view -H -m2 -M2 -v snps $W1943 2>/dev/null | wc -l
echo "new BCF records:"; bcftools view -H $OUTGRP/W1943.bcf 2>/dev/null | wc -l
bcftools query -l $OUTGRP/W1943.bcf 2>/dev/null
}

#RENAME OUTGROUP CHR AND SPLIT BY CHR
function prep_outgrp {

OG=$1        # e.g. $OUTGRP/W1715.bcf  (already biallelic-SNP filtered + reheadered)
NAME=$2      # W1715 or W1943
OUTDIR=$3

source $HOME/miniforge3/bin/activate est-sfs

cat > rename.chr.toNum.txt <<'EOF'
chr01	1
chr02	2
chr03	3
chr04	4
chr05	5
chr06	6
chr07	7
chr08	8
chr09	9
chr10	10
chr11	11
chr12	12
EOF

# 1. rename chr01->1 ... chr12->12, keep only chr1-12, write bgzipped VCF
bcftools annotate --rename-chrs rename.chr.toNum.txt $OG -x INFO,'^FORMAT/GT' -Ou \
  | bcftools view -t 1,2,3,4,5,6,7,8,9,10,11,12 -Oz -o $OUTDIR/${NAME}.renamed.vcf.gz
tabix -p vcf $OUTDIR/${NAME}.renamed.vcf.gz

# 2. split by chromosome
for i in $(seq 1); do
    bcftools view $OUTDIR/${NAME}.renamed.vcf.gz -r $i \
        -Oz -o $OUTDIR/${NAME}.chr${i}.vcf.gz
    tabix -p vcf $OUTDIR/${NAME}.chr${i}.vcf.gz
done

# 3. verify
echo "$NAME chroms:"; tabix -l $OUTDIR/${NAME}.renamed.vcf.gz
}

#CHECK IF THE REFERENCES FOR FOCAL AND OUTGROUPS ARE THE SAME
function check_vcf {

DIR1=$1 
DIR2=$2         # dir holding the per-chr vcf.gz files
FOCALPREFIX=$3  # e.g. genome1.focal   -> $DIR/genome1.focal.chr${i}.vcf.gz
OG1=$4          # e.g. W1943           -> $DIR/W1943.chr${i}.vcf.gz
OG2=$5          # e.g. W1715           -> $DIR/W1715.chr${i}.vcf.gz

source $HOME/miniforge3/bin/activate est-sfs

for i in $(seq 1 2); do

    FOCAL=$DIR1/${FOCALPREFIX}.chr${i}.vcf.gz
    OUT1=$DIR2/${OG1}.chr${i}.vcf.gz
    OUT2=$DIR2/${OG2}.chr${i}.vcf.gz

    echo "====================checking vcf chr${i} ===================="

    # --- 0. files exist? ---
    for f in $FOCAL $OUT1 $OUT2; do
        if [ ! -s "$f" ]; then
            echo "  MISSING or empty: $f"; continue 2
        fi
        # make sure indexed (needed for -r region queries)
        tabix -p vcf "$f" 2>/dev/null
    done

    # --- 1. chromosome name actually present in each file ---
    FC=$(tabix -l $FOCAL 2>/dev/null | tr '\n' ',')
    O1C=$(tabix -l $OUT1 2>/dev/null | tr '\n' ',')
    O2C=$(tabix -l $OUT2 2>/dev/null | tr '\n' ',')
    echo "  chrom in focal: $FC | $OG1: $O1C | $OG2: $O2C"
    if [ "$FC" != "${i}," ] || [ "$O1C" != "${i}," ] || [ "$O2C" != "${i}," ]; then
        echo "  WARNING: chromosome naming not uniformly '${i}' -- check before merge!"
    fi

    # --- 2. contig lengths match (same reference build) ---
    bcftools view -h $FOCAL 2>/dev/null | grep "^##contig" | grep -E "ID=${i}[,>]" | sort > contigs.focal.$i.txt
    bcftools view -h $OUT1  2>/dev/null | grep "^##contig" | grep -E "ID=${i}[,>]" | sort > contigs.out1.$i.txt
    bcftools view -h $OUT2  2>/dev/null | grep "^##contig" | grep -E "ID=${i}[,>]" | sort > contigs.out2.$i.txt
    echo "  -- contig length focal vs $OG1 --"; diff contigs.focal.$i.txt contigs.out1.$i.txt && echo "     identical"
    echo "  -- contig length focal vs $OG2 --"; diff contigs.focal.$i.txt contigs.out2.$i.txt && echo "     identical"

    # --- 3. REF-allele consistency at shared positions ---
    bcftools query -f '%CHROM\t%POS\t%REF\n' $FOCAL 2>/dev/null | awk '{print $1"_"$2"\t"$3}' | sort > ref.focal.$i.txt
    bcftools query -f '%CHROM\t%POS\t%REF\n' $OUT1  2>/dev/null | awk '{print $1"_"$2"\t"$3}' | sort > ref.out1.$i.txt
    bcftools query -f '%CHROM\t%POS\t%REF\n' $OUT2  2>/dev/null | awk '{print $1"_"$2"\t"$3}' | sort > ref.out2.$i.txt

    echo -n "  $OG1: "; join -j 1 ref.focal.$i.txt ref.out1.$i.txt \
        | awk '{tot++; if($2!=$3) mm++} END{printf "shared=%d  REF_mismatch=%d  (%.4f%%)\n", tot, mm, (tot? 100*mm/tot:0)}'
    echo -n "  $OG2: "; join -j 1 ref.focal.$i.txt ref.out2.$i.txt \
        | awk '{tot++; if($2!=$3) mm++} END{printf "shared=%d  REF_mismatch=%d  (%.4f%%)\n", tot, mm, (tot? 100*mm/tot:0)}'

    # cleanup temp files for this chr
    rm -f contigs.focal.$i.txt contigs.out1.$i.txt contigs.out2.$i.txt \
          ref.focal.$i.txt ref.out1.$i.txt ref.out2.$i.txt
done
}

#MERGE FOCAL AND OUTGROUP VCF
#MERGE FOCAL AND OUTGROUPS, PER CHROMOSOME
function merge_vcf {

DIR1=$1  
DIR2=$2        # dir holding the per-chr vcf.gz files
FOCALPREFIX=$3  # e.g. genome1.focal  -> $DIR/genome1.focal.chr${i}.vcf.gz
OG1=$4          # e.g. W1943          -> $DIR/W1943.chr${i}.vcf.gz
OG2=$5          # e.g. W1715          -> $DIR/W1715.chr${i}.vcf.gz
OUTDIR=$6       # where to write merged files

source $HOME/miniforge3/bin/activate est-sfs
echo "USING ENV" $CONDA_DEFAULT_ENV

for i in $(seq 1); do

    FOCAL=$DIR1/${FOCALPREFIX}.chr${i}.vcf.gz
    OUT1=$DIR2/${OG1}.chr${i}.vcf.gz
    OUT2=$DIR2/${OG2}.chr${i}.vcf.gz

    echo "==================== merging chr${i} ===================="

    # --- 0. sanity: inputs exist and are indexed (merge needs indexes) ---
    for f in $FOCAL $OUT1 $OUT2; do
        if [ ! -s "$f" ]; then
            echo "  MISSING or empty: $f -- skipping chr${i}"; continue 2
        fi
        [ -s "${f}.tbi" ] || tabix -p vcf "$f"
    done

    # --- 1. merge the three files ---
    bcftools merge --threads 4 \
        $FOCAL \
        $OUT1 \
        $OUT2 \
        -Oz -o $OUTDIR/merged.tmp.chr${i}.vcf.gz 
    #treat missing data in og as ref: not ./. but 0/0, key for og merging
    bcftools +setGT $OUTDIR/merged.tmp.chr${i}.vcf.gz -Oz \
        -o $OUTDIR/merged.chr${i}.vcf.gz \
        -- --target-gt . \
        --new-gt 0 \
        --samples W1943,W1715

    tabix -p vcf $OUTDIR/merged.chr${i}.vcf.gz

    # --- 3. verify ---
    echo "  samples:"; bcftools query -l $OUTDIR/merged.chr${i}.vcf.gz 2>/dev/null | tr '\n' ' '; echo
    echo "  chrom:"; tabix -l $OUTDIR/merged.chr${i}.vcf.gz
    echo "  SNPs (final):";       bcftools view -H $OUTDIR/merged.chr${i}.vcf.gz 2>/dev/null | wc -l

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
bcftools query -l $INFILE |grep -v -E '^(W1943|W1715)$' > focal.txt 

python3 prep_estsfs.py \
    --infile $INFILE \
    --focalSampleList focal.txt \
    --outgroup1SampleList og1.txt \
    --outgroup2SampleList og2.txt \
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
awk -v OFS=',' '
    # Load outgroup file first (NR==FNR reads first file)
    NR==FNR { 
        og[$1] = $2 #   Store a hash: og[chrom:pos] = base
        next
    }
    # Process reference file
    {
        key = $1":"$2          # Build a lookup key from ref: "chrom:pos"e.g. "1:1118"
        base = og[key]         # look up outgroup base
        if (base == "A")      print 1,0,0,0
        else if (base == "C") print 0,1,0,0
        else if (base == "G") print 0,0,1,0
        else if (base == "T") print 0,0,0,1
        else                  print 0,0,0,0   # missing in outgroup
    }
' $OG_FILE $REF_FILE > $OUTDIR/allele_count_$OGNAME.txt

echo "Done. Output written to $OUTDIR/allele_count_$OGNAME.txt"
echo "Total positions: $(wc -l < $OUTDIR/allele_count_$OGNAME.txt)"

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

#WGA WITH CACTUS FOR OUTLIER ALLELE EXTRACTION
#run_repeatmasker $GLAMPTL
#run_cactus $IRGSP $BARTHII $GLAMPTL $RESULT1 Obrth Oglmptl Ostv
#get_vcf_coordinates $VCF_ALL $RESULT1
#get_hal_coordinates $RESULT1/rice_wga.hal $RESULT1/ref_snps.bed Ostv Obrth $RESULT1
#get_hal_coordinates $RESULT1/rice_wga.hal $RESULT1/ref_snps.bed Ostv Oglmptl $RESULT1  
#get_OG_bases $RESULT1/snps_Obrth.bed $BARTHII Obrth $RESULT1
#get_OG_bases $RESULT1/snps_Oglmptl.bed $GLAMPTL Oglmptl $RESULT1
#sfs_OG $RESULT1/snps_IRGSP.tsv $RESULT1/snps_Obrth.tsv Obrth $RESULT1
sfs_OG $RESULT1/snps_IRGSP.tsv $RESULT1/snps_Oglmptl.tsv Oglmptl $RESULT1


#CONVERT PLINK TO VCF
#plink_to_vcf $SAMPLE $RESULT1 genome1.biallelic.base.genomewide.SNPs.withID.PlinkFormat genome1.focal

#MAKE SURE ONLY INCLUDE BIALLILIC SNPs, rename chr and extract chr 1-12
#biallelic $RESULT1/genome1.focal.vcf.gz $RESULT1/genome1.focal.biallelic.vcf.gz




#rename_outgrp #DONE
#prep_outgrp $OUTGRP/W1943.bcf W1943 $OUTGRP
#prep_outgrp $OUTGRP/W1715.bcf W1715 $OUTGRP

#CHECK IF THE REF IS THE SAME AS THE OUTGROUP VCF
#check_vcf $RESULT1 $OUTGRP genome1.focal W1943 W1715

#MERGE FOCAL AND OUTGROUP
#merge_vcf $RESULT1 $OUTGRP genome1.focal W1943 W1715 $RESULT1

#GENERATE EST-SFS INPUT FILES
#count_allele $RESULT1/merged.chr1.vcf.gz $RESULT1 chr1

#RUN EST-SFS
#run_est_sfs $RESULT1/estsfs_data.chr1.tsv $R6 $RESULT1 chr1
#run_est_sfs chr1_first100.tsv $R6 $RESULT1 chr1
}


main