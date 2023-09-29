#!/bin/bash
#
# prepare all the BCF files required for phasing -  minimal so as to be used within the snakemake pipeline
# should work given a list of trios and the QCed genotypes. It expects one file for WES+Array, and another for all rare. 
#

module load common-apps/bcftools/1.16

### input ###
chr=$1
tag=$2
work_dir=$3
pedigree="$work_dir/$tag.pedigree"
input_rare="FIXTHIS/chr${chr}_QCed.vcf.gz"
input_merged="FIXTHIS/chr${chr}.merged_WES_GSA.vcf.gz"

### output ###
out_bcf_common_full="${work_dir}/sandbox/${tag}.notrios.merged_WES_GSA.chr${chr}.bcf"
out_bcf_common_trios="${work_dir}/sandbox/${tag}.trios.merged_WES_GSA.chr${chr}.bcf"
out_bcf_rare_full="${work_dir}/sandbox/${tag}.notrios.rare_prepared.chr${chr}.bcf"
out_bcf_rare_trios="${work_dir}/sandbox/${tag}.trios.rare_prepared.chr${chr}.bcf"

if [ ! -f $work_dir/$tag.trios.parents ]; then 
    echo "Generate a list of $tag and the corresponding parents"
    # NOTE: This is now performed with python according to mendel-errors
    cut -f2,3 $pedigree | tr "\t" "\n" |  sort | uniq > $work_dir/$tag.trios.parents 
    cat $pedigree | tr "\t" "\n" | sort  | uniq > $work_dir/$tag.trios.samples
fi
echo "Samples to keep for trio-phasing:" $(wc -l $work_dir/$tag.samples | awk '{print $1}')
echo "Generating new BCFs for common variants..."

SECONDS=0
# keep only the samples in the trio cohort
bcftools view $input_merged -S $work_dir/$tag.trios.samples -Ob -o $out_bcf_common_trios
bcftools index $out_bcf_common_trios
echo "Done with preprocessing common-trios, duration: ${SECONDS}."
SECONDS=0
# keep all but the parents in the trio cohort
bcftools view $input_merged -S ^$work_dir/$tag.trios.parents -Ob -o $out_bcf_common_full
bcftools index $out_bcf_common_full
echo "Done with preprocessing common-full, duration: ${SECONDS}."

if [ ! -f $work_dir/$tag.WES_final.samples ]; then
    bcftools query $out_bcf_common_full --list samples > $work_dir/$tag.WES_final.samples
fi
echo "Samples to keep for downstream analyses:" $(wc -l $work_dir/$tag.WES_final.samples | awk '{print $1}')

echo -e "\n################\n"
echo "Generating new BCFs for rare variants..."
# here we have to first change some sample IDs, then select samples, then exclude variants of higher freq
# we had to change a few sample IDs in which case we'd need 
# bcftools reheader $input_rare --samples $work_dir/samples.update_ids_wes.txt | bcftools view -S $work_dir/$tag.samples -Ou | bcftools filter --exclude 'MAF>0.001' -Ob -o $out_bcf_rare_full
SECONDS=0
bcftools view $input_rare -S $work_dir/$tag.trios.samples -Ou | bcftools filter --exclude 'MAF>0.001' -Ob -o $out_bcf_rare_trios
bcftools index $out_bcf_rare_trios
echo "Done with preprocessing rare-trios, duration: ${SECONDS}."
SECONDS=0
bcftools view $input_rare -S $work_dir/$tag.WES_final.samples -Ou | bcftools filter --exclude 'MAF>0.001' -Ob -o $out_bcf_rare_full
bcftools index $out_bcf_rare_full
echo "Done with preprocessing rare-trios, duration: ${SECONDS}."

echo "All done!"
