#!/bin/bash

final_prefix=$1
chr_vcf_prefix=$2

files_to_concat=$1.files.txt
rm -f $files_to_concat
for C in {1..22}; do echo $chr_vcf_prefix.chr$C.vcf >> $files_to_concat; done
bcftools concat -f $files_to_concat -Ob -o $final_prefix.bcf

plink2 --export bgen-1.2 bits=8 --bcf $final_prefix --out $final_prefix
# almost done, fix the .sample file
head -2 $final_prefix.sample > $final_prefix.sample1
awk 'NR>2{print 1,$2,0,0}' $final_prefix.sample >> $final_prefix.sample1
mv $final_prefix.sample1 $final_prefix.sample