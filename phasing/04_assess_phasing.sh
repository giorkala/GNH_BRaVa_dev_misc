#!/bin/bash

# GK - Jul 28th, 2023
# based on https://odelaneau.github.io/shapeit5/docs/documentation/switch/
# note, all files should be BCF for higher efficiency
# see README.md for more details

if [[ -z $1 ]]; then
    chr=$LSB_JOBINDEX
else
    chr=$1
fi

SHAPEIT_switch='/software/team281/bin/shapeit5/switch_static'
threads=5

tag=100trios
work_dir="/FIXTHIS/phasing"
mkdir -p $work_dir/phasing_assessment/

###############################################
echo "Work for chrom-$chr and common variants:"

phased_true="${work_dir}/phased_genotypes_common/GNH_39k.chr${chr}.phased_common.$tag.bcf"
phased_estd="${work_dir}/phased_genotypes_common/GNH_39k.chr${chr}.phased_common.no$tag.bcf"
out_prefix="$work_dir/phasing_assessment/assess.common.chr$chr"

# $SHAPEIT_switch --validation $phased_true --estimation $phased_estd --region chr$chr --output $out_prefix --thread $threads

zcat $out_prefix.sample.switch.txt.gz | awk 'BEGIN { e=0; t=0; } { e+=$2; t+=$3; } END { print "SER =", e*100/t; }'

#############################################
echo -e "\n################\n"
echo "Work for chrom-$chr and rare variants:"

phased_true="${work_dir}/phased_genotypes_rare/GNH_39k.chr${chr}.phased_rare.$tag.bcf"
phased_estd="${work_dir}/phased_genotypes_rare/GNH_39k.chr${chr}.phased_rare.no$tag.bcf"
out_prefix="$work_dir/phasing_assessment/assess.rare.chr$chr"

for PP in 0.50 0.90 0.95 0.99; do 
    out_prefix="$work_dir/phasing_assessment/assess.rare.pp$PP.chr$chr"
    $SHAPEIT_switch --validation $phased_true --estimation $phased_estd --min-pp $PP \
    --singleton \
    --region chr$chr --output $out_prefix --thread $threads --log $out_prefix.log
done

# $SHAPEIT_switch --validation $phased_true --estimation $phased_estd --region chr$chr --output $out_prefix --thread $threads

zcat $out_prefix.sample.switch.txt.gz | awk 'BEGIN { e=0; t=0; } { e+=$2; t+=$3; } END { print "SER =", e*100/t; }'


# run as:
# bsub -J assess_phasing[20-22] -o logs/assess_phasing.%I.%J -q normal -R"select[mem>4000] rusage[mem=4000]" -M4000 bash 04_assess_phasing.sh
