#!/usr/bin/bash

# This script will survey each gene and prepare one VCF with the genotypes for biallelic samples.
# The script should be called by the Snakefile in parallel for each gene.

## INPUT ##
BCF=$1 # original BCF with phased genotypes
annot=$2 # annotation file
GENE=$3 # gene name
CHR=$4 # chromosome
gene_final=$5
# work_dir=$6 # working directory
gene_tmp_prefix=$6/sandbox/tmp.$GENE # prefix for temp files

# required files assumed at hand:
samples_all=$6/samples.bcf.txt
samples_biallelic=$6/sandbox/samples.$GENE.biallelic
samples_other=$6/sandbox/samples.$GENE.other

# Check if the input argument matches the pattern "chr" followed by a number from 1 to 22
if ! [[ $CHR =~ ^chr([1-9]|1[0-9]|2[0-2])$ ]]; then
    echo "Warning: invalid chrom index! Switching to chr$CHR"
    CHR="chr$CHR"
fi

# generate a list of suitable variants
variants="$gene_tmp_prefix.variants"
rm -rf $variants
# note: the following grep-ing is reduntant as we've already pre-selected variants
for consq in pLoF damaging_missense_or_protein_altering other_missense_or_protein_altering; do # without 
    grep -w $GENE $annot | grep $consq | awk '{print $1}' >> $variants
done
n_v=$(wc -l $variants | awk '{print $1}')
# pull out the region with observed variants for the current gene
R1=$(grep $GENE $annot | cut -d: -f2 | head -1)
R2=$(grep $GENE $annot | cut -d: -f2 | tail -n1)
echo "Variants found for $CHR-$GENE: $n_v; Region: $R1 - $R2"

# make ONE VCF file will all samples but gene-specific variants; this will speed up the process
# NOTE: you might wanna add SNP IDs in this step using `bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT'`
bcft_cmd="bcftools view -r $CHR:$R1-$R2 $BCF --threads 2 -Oz -o $gene_tmp_prefix.all.vcf.gz"
eval $bcft_cmd
echo "Gene-specific VCF is generated, now subsampling to biallelic individuals..."
BCF="$gene_tmp_prefix.all.vcf.gz"

# now make two VCF files, one per set, then merge them
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
