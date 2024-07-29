#!/usr/bin/bash

# This script will survey each gene and prepare one VCF with the genotypes for biallelic samples.
# The script should be called by the Snakefile in parallel for each gene.

## INPUT ##
BCF=$1 # original BCF with phased genotypes
annot=$2 # annotation file
GENE=$3 # gene name
CHR=$4 # chromosome
gene_final=$5
gene_tmp_prefix=$6/sandbox/tmp.$GENE # prefix for temp files

# requirec files assumed at hand:
samples_all=$6/samples.bcf.txt
samples_biallelic=$6/sandbox/samples.$GENE.biallelic
samples_other=$6/sandbox/samples.$GENE.other

# generate a list of suitable variants
variants="$gene_tmp_prefix.variants"
rm -rf $variants
# note: the following grep-ing is reduntant as we've already pre-selected variants
for consq in pLoF damaging_missense_or_protein_altering other_missense_or_protein_altering; do # without 
    grep -w $GENE $annot | grep $consq | cut -f1 >> $variants
done
n_v=$(wc -l $variants | awk '{print $1}')
echo "Variants found for $GENE: $n_v"

# make two VCF files, one per set, and merge them
# select samples and variants from a VCF file - I'm using `eval` as we need to pass a variable in the filtering argument
bcft_cmd="bcftools view $BCF -i'ID=@$variants' -S $samples_biallelic -Oz -o $gene_tmp_prefix.biallelic.vcf.gz"
eval $bcft_cmd
bcftools index $gene_tmp_prefix.biallelic.vcf.gz 
bcft_cmd="bcftools view $BCF -i'ID=@$variants' -S $samples_other -Oz -o $gene_tmp_prefix.other.vcf.gz"
eval $bcft_cmd
# now turn all genotypes to 0/0 for the latter
bcftools filter $gene_tmp_prefix.other.vcf.gz -i'GT="1/0 | 0/1"' --set-GT 0 -Oz -o $gene_tmp_prefix.other_zeroed.vcf.gz
# echo "check"
bcftools index $gene_tmp_prefix.other_zeroed.vcf.gz
echo "Merging the two files and finalising..."
# finally, merge the two files and reorder the samples to enable concatenation later -- | bcftools reheader -s (this didnt work...)
bcftools merge $gene_tmp_prefix.biallelic.vcf.gz $gene_tmp_prefix.other_zeroed.vcf.gz -Ou | bcftools view -c 1 -S $samples_all -Ob -o $gene_final

# all done! delete temp files
rm $gene_tmp_prefix.* #$tmp_dir/samples.$GENE.* 
echo "All done for $GENE. Duration = ${SECONDS} (secs)."