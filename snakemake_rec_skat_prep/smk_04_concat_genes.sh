#!/bin/bash

# create a chrom-wide VCF using all gene-level files
# to do so, we need to identify all available genes (=having a VCF) for a given chrom
# then add them to a single VCF file with `bcftools concat ...`
# then remove the PP field and turn any phased genotype to unphased

## INPUT ##
CHR=$1
genes_list=$2
out_prefix=$3
# this should be "$work_dir/vcf_files_temp/$TAG.chr$CHR"
out_final="$out_prefix.vcf"

# temp files needed
files_to_concat="$out_prefix.files_used"
files_missed="$out_prefix.files_missed"
rm -f $files_to_concat $files_missed

# first, check which files are available
cat $genes_list | while read chrom GENE; do
    candidate="$out_prefix.$GENE.vcf"
    if [ -f $candidate ]; then
        echo $candidate >> $files_to_concat
    else
        echo $candidate >> $files_missed
    fi
done
echo "Genes available for chr$CHR: $(wc -l $files_to_concat | awk '{print $1}')"

# now merge all gene-VCF files to one
bcftools concat -f $files_to_concat -Ou | bcftools sort -m 2G -Ou | bcftools norm --rm-dup all -Ob -o $out_prefix.bcf
# make a new file (NOTE: might wanna simplify the header)
bcftools view -h $out_prefix.bcf > $out_final
# bcftools view -h $out_prefix.bcf | tail -n17 >> $out_final
# then drop the PP and turn any phased genotype to unphased
bcftools annotate -x FORMAT/PP $out_prefix.bcf -Ou | bcftools view -H | sed 's/|/\//g' >> $out_final

echo "Compression and indexing..."
bgzip -c $out_final
bcftools index $out_final.gz

rm $out_prefix.bcf
