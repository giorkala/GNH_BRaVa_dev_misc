#!/usr/bin/bash

## INPUT ##
BCF=$1 # BCF with phased genotypes {input.input_bcf} 
annot=$2 # annotation file {input.annotation} Regenie-like [chr:pos:ref:alt gene consequence]
max_AF=$3
WORK_DIR=$4
TAG=$5
output=$6
# there's an additional parameter: min_PP=0.90 at line 37 (the awk in the pipe)

tmp_dir="$WORK_DIR/sandbox"
genotypes="$tmp_dir/alt_gen.$TAG.PP90af$max_AF.txt.gz"

call_chets="/software/team281/gk18/call_chets/call_chets"

echo "Will create a file with biallelic carriers for $TAG."

if [ ! -f $genotypes ]; then 
    echo "Generating a file with all genotypes. This can be slow..."
    # select all heterozygous genotypes with PP>0.90 or missing (which implies a common variant)
    bcftools view $BCF --max-af $max_AF | bcftools query -i '(GT="1|1" | GT="0|1" | GT="1|0")' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT %PP\n]' | awk '(length($2)<50)' | awk '($4=="." || $4>0.90)' | gzip > $genotypes
    # what is wrong here? bcftools view $BCF --max-af 0.05 -Ou | bcftools query -f '[%CHROM:%POS %SAMPLE %GT\n]' -i'GT="alt"' |
fi

# call the C++ program to identify carriers and filter to homozygotes and comp-hets
$call_chets -g $genotypes -m $annot | grep 'chet\|hom' > $output

# end of step1