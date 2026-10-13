#!/bin/bash
function set_up_env_ancify {

CONDA_SUBDIR=osx-64 conda create -y -n ancify-env python=3.9
#activate env
source $HOME/miniforge3/bin/activate ancify-env
echo $CONDA_DEFAULT_ENV #double check the env is correct

#permanently locks that particular Conda environment to only accept osx-64 (Intel) packages
conda config --env --set subdir osx-64
#install
conda install -c bioconda lastz -n ancify-env -y
conda install -c bioconda ucsc-axtchain -n ancify-env -y
conda install -c bioconda ucsc-chainmergesort -n ancify-env -y
conda install -c bioconda ucsc-chainprenet -n ancify-env -y
conda install -c bioconda ucsc-chainnet -n ancify-env -y
conda install -c bioconda ucsc-nettoaxt -n ancify-env -y
conda install -c bioconda ucsc-netsyntenic -n ancify-env -y
conda install -c bioconda ucsc-netchainsubset -n ancify-env -y
conda install -c bioconda ucsc-fatotwobit ucsc-twobitinfo -n ancify-env -y
conda install -c bioconda ucsc-chainstitchid ucsc-maftoaxt -n ancify-env -y
conda install -c bioconda bx-python -n ancify-env -y
}

#CONVERT FASTA TO 2BIT

function convert_fa_to_2bit {

    INFA=$1
    NAME=$2
    OUTDIR=$3

source $HOME/miniforge3/bin/activate ancify-env
echo "USING ENV" $CONDA_DEFAULT_ENV

echo "CONVERTING FASTA TO 2bit"
echo "INFILE: ${INFA}"
echo "OUTFILE: ${OUTDIR}/${NAME}.2bit"

    fatotwobit $INFA ${OUTDIR}/${NAME}.2bit
    twoBitInfo ${OUTDIR}/${NAME}.2bit ${OUTDIR}/${NAME}.sizes
}

#use the available hal file:
function hal_2_maf {

HAL=$1
NAME=$2
OUTDIR=$3

source $HOME/miniforge3/bin/activate cactus-env
echo "USING ENV" $CONDA_DEFAULT_ENV

echo "convert hal to maf: ${NAME}"

hal2maf $HAL $OUTDIR/${NAME}.maf \
    --refGenome Ostv \
    --targetGenomes ${NAME} \
    --noAncestors \
    --maxRefGap 100
}

function MAF_2_AXT {

MAF=$1
T_NAME=$2
Q_NAME=$3
OUTDIR=$4

source $HOME/miniforge3/bin/activate ancify-env
echo "USING ENV" $CONDA_DEFAULT_ENV

echo "convert MAF to AXT: ${QRY_NAME}"

#maftoaxt $MAF first $Q_NAME $OUTDIR/${Q_NAME}.axt
maf_to_axt.py $T_NAME $Q_NAME <$MAF> $OUTDIR/${Q_NAME}.axt
}

#
function AXT_2_Net_AXT {

Q=$1
T=$2
OUTDIR=$3

source $HOME/miniforge3/bin/activate ancify-env
echo "USING ENV" $CONDA_DEFAULT_ENV

echo "convert AXT 2 Net-AXT: ${NAME}"
# 1. Convert AXT to Chain
axtChain -linearGap=medium $OUTDIR/${Q}.axt $OUTDIR/${T}.2bit $OUTDIR/${Q}.2bit $OUTDIR/${Q}.chain

# 2. Merge and sort chains
chainMergeSort $OUTDIR/${Q}.chain > $OUTDIR/${Q}.all.chain

# 3. PreNet (removes chains with low coverage)
chainPreNet $OUTDIR/${Q}.all.chain $OUTDIR/${T}.sizes $OUTDIR/${Q}.sizes $OUTDIR/${Q}.preNet.chain

# 4. Create Net: usage: chainNet in.chain target.sizes query.sizes target.net query.net
#use output over the target genome
chainNet $OUTDIR/${Q}.preNet.chain $OUTDIR/${T}.sizes $OUTDIR/${Q}.sizes $OUTDIR/${T}-${Q}.net $OUTDIR/${Q}.net

# 5. Filter to Syntenic Nets
#netSyntenic $OUTDIR/${T}-${Q}.net $OUTDIR/${T}-${Q}.syn.net #optional: syntex annotation
netChainSubset $OUTDIR/${T}-${Q}.net $OUTDIR/${Q}.preNet.chain $OUTDIR/${T}-${Q}.net.chain

# 6. Stitch duplicate chain IDs together
#echo "STITCH DUPLICATE CHAIN IDs TOGETHER"
chainStitchId $OUTDIR/${T}-${Q}.net.chain $OUTDIR/${T}-${Q}.net.stitched.chain

# 7. Convert Net to AXT:  netToAxt in.net in.chain target.2bit query.2bit out.axt
netToAxt $OUTDIR/${T}-${Q}.net $OUTDIR/${T}-${Q}.net.stitched.chain  $OUTDIR/${T}.2bit $OUTDIR/${Q}.2bit $OUTDIR/${T}-${Q}.net.axt

# 8. Compress the final AXT
gzip $OUTDIR/${T}-${Q}.net.axt

rm -f $OUTDIR/${Q}.chain $OUTDIR/${Q}.all.chain $OUTDIR/${Q}.preNet.chain $OUTDIR/${T}-${Q}.net $OUTDIR/${Q}.net $OUTDIR/${T}-${Q}.net.chain ${T}-${Q}.net.stitched.chain $OUTDIR/${T}-${Q}.syn.net
}
#---------------------execute------------------------

function main {  

OUTGRP=/Users/zchen13/data/outgroups
SAMPLE=/Users/zchen13/data/oryza_19K_RGP
IRGSP=/Users/zchen13/data/reference/IRGSP-1.0_chr.renamed.fa
BARTHII=/Users/zchen13/data/reference/O_barthii.chr.fa
GLAMPTL=/Users/zchen13/data/reference/O_glamaepatula.chr.fa
MRDONLS=/Users/zchen13/data/reference/O_meridionalis.chr.fa
WORKDIR=/Users/zchen13/SCRATCH/ancestral_allele
#
RESULT1=$WORKDIR/01_est-sfs
RESULT2=$WORKDIR/02_est-sfs
RESULT3=$WORKDIR/03_ancify #AA INFERENCE USING ANCIFY


#OUTPUT SHORTCUTS
VCF_ALL=$RESULT1/genome1.focal.chr1-12.vcf.gz

#est-sfs config files
JC=/Users/zchen13/scripts/ancestral-allele/config-JC.txt
KIMURA=/Users/zchen13/scripts/ancestral-alleleconfig-kimura.txt
R6=/Users/zchen13/scripts/ancestral-allele/config-rate6.txt

#set_up_env_ancify

function setup_dir {

mkdir -p $WORKDIR
#mkdir -p $RESULT1
#mkdir -p $RESULT2
mkdir -p $RESULT3
#mkdir -p $RESULT5
#mkdir -p $RESULT5

}
setup_dir

#CONVERT FA --> 2BIT
#convert_fa_to_2bit $IRGSP IRGSP $RESULT3 
#convert_fa_to_2bit $BARTHII Obrth $RESULT3 
#convert_fa_to_2bit $GLAMPTL Oglmptl $RESULT3 
#convert_fa_to_2bit $MRDONLS Omrdonls $RESULT3 

#HAL 2 MAF
#hal_2_maf $RESULT1/rice_wga.hal Obrth $RESULT3
#hal_2_maf $RESULT1/rice_wga.hal Oglmptl $RESULT3
#hal_2_maf $RESULT2/rice_wga.hal Omrdonls $RESULT3

#MAF 2 AXT
#MAF_2_AXT $RESULT3/Obrth.maf    Ostv Obrth    $RESULT3
#MAF_2_AXT $RESULT3/Oglmptl.maf  Ostv Oglmptl  $RESULT3
MAF_2_AXT $RESULT3/Omrdonls.maf Ostv Omrdonls $RESULT3

#AXT 2 NET AXT
#AXT_2_Net_AXT Obrth IRGSP $RESULT3
#AXT_2_Net_AXT Oglmptl IRGSP $RESULT3
AXT_2_Net_AXT Omrdonls IRGSP $RESULT3

#RUN ANCIFY333
}


main