#!/bin/bash
#
# Merge WES data with array from the GSA cohort
# this is likely to fail as a lot of I/O is required, so please modify accordingly
# particularly for step-1, as we had to consider several changes for sample IDs which might not be needed.
# 
# based on https://github.com/frhl/wes_ko_ukbb/blob/main/scripts/phasing/prefilter/06_wes_union_calls.sh
# see README.md for more details

module load common-apps/bcftools/1.16
plink=/software/team281/bin/plink # we could provide cpu/ram limits for this

### input ###
chr=$1
tag=$2
work_dir=$3

path_array=FIXTHIS/GSA_QC
prefix_array=$path_gsa/GNH_GSAv3EAMD_DupExcl_autosome_atcgSNP
# note: we assumme one BED for all autosomes, but that could be easily changed
path_wes=FIXTHIS/filtered_vcfs
prefix_wes=${path_wes}/chr${chr}_hard_filters.vcf.gz

### output ###
wes_new_prefix=$work_dir/sandbox/tmp.chr$chr.wes_filtered
gsa_new_prefix=$work_dir/sandbox/tmp.chr$chr.gsa_filtered
final_out=$work_dir/sandbox/chr$chr.merged_WES_GSA.vcf.gz


## 1. make lists for the intersection of samples ##
SECONDS=0
# for that, I first focus on WES who are in GSA, and then refine the sample according to who's left
# this can be simplified
if [ ! -f $work_dir/samples.WES_in_GSA.fam ]; then
    # create an index of samples in the WES cohort
    bcftools query $wes_prefix --list samples > $work_dir/samples.WES_all.txt
    # find the intersection of those with the GSA cohort - the link file contains more
    join -j 1 $work_dir/samples.WES_all.txt <(awk 'NR>1{if($2>0){print $0} }' $work_dir/link_OrageneID_all-WES_GSA.txt | sort) > $work_dir/samples.common_WES_GSA.link 
    # for t in 2 3 4; do sed -i -E "s/\_$t//g" $work_dir/samples.common_WES_GSA.link; done
    awk '{print $1}' $work_dir/samples.common_WES_GSA.link > $work_dir/samples.WES_in_GSA.txt
    # also create a new linkage file according to PLINK's required format for updating IDs
    join -1 2 -2 1 <(cat samples.common_WES_GSA.link | sort -t' ' -k2 ) <(cat $gsa_prefix.fam | sort) | awk '{print $1,$3,$2,$2}' > $work_dir/samples.update_ids_gsa.txt
    sed -i -e 's/\_2//g' -e 's/\_3//g' -e 's/\_4//g' $work_dir/samples.update_ids_gsa.txt
    awk '{print $1,$1}' $work_dir/samples.common_WES_GSA.link | sed -e 's/\_2//g' -e 's/\_3//g' -e 's/\_4//g' > $work_dir/samples.WES_in_GSA.fam 
fi

## 2. make a new VCF focusing on the new set of samples and MAF>0.1% ##
if [ ! -f $wes_new_prefix.vcf.gz ]; then
    echo -e "\nCalling BCFtools to prepare the WES input. This is will be slow!"
    bcftools view $wes_prefix -S $work_dir/samples.WES_in_GSA.txt -Ou | bcftools reheader --samples $work_dir/samples.update_ids_wes.txt | bcftools filter --exclude 'MAF<0.001' -Oz -o $wes_new_prefix.vcf.gz
fi

## 3. prepare the GSA genotypes ##
# focusing on the new set of samples (on the original index)
# subselect individuals, change from the original index to the VCF one, replace "22" with "chr22", and make the vcf in one go
echo -e "\nCalling PLINK to prepare the GSA input."
$plink --bfile $gsa_prefix --chr $chr \
    --maf 0.001 \
    --update-ids $work_dir/samples.update_ids_gsa.txt \
    --keep $work_dir/samples.WES_in_GSA.fam \
    --real-ref-alleles \
    --recode vcf-fid bgz \
    --output-chr chr26 \
    --out $gsa_new_prefix

## 4. make a list of overlapping variants to exclude from GSA ##
bcftools view -HG $wes_new_prefix.vcf.gz | awk '{ print $1":"$2}' > $wes_new_prefix.snpid
bcftools view -HG $gsa_new_prefix.vcf.gz | awk '{ print $1":"$2}' > $gsa_new_prefix.snpid
comm -23 <( sort $gsa_new_prefix.snpid) <( sort $wes_new_prefix.snpid) > $gsa_new_prefix.snpid_unique
sed -i 's/:/\t/g' $gsa_new_prefix.snpid_unique
a=$(wc -l $gsa_new_prefix.snpid | awk '{print $1}')
b=$(wc -l $gsa_new_prefix.snpid_unique | awk '{print $1}')
echo "Overlapping variants to remove:" $((a-b))

## 5. sort individuals in GSA's VCF according to the former ##
echo -e "\nProceeding with BCFtools and final merging."
bcftools query $wes_new_prefix.vcf.gz --list samples > $work_dir/samples.WES.final
bcftools index $gsa_new_prefix.vcf.gz
bcftools index $wes_new_prefix.vcf.gz
bcftools view $gsa_new_prefix.vcf.gz -S $work_dir/samples.WES.final -R $gsa_new_prefix.snpid_unique -Oz -o $gsa_new_prefix.sorted.vcf.gz
bcftools index $gsa_new_prefix.sorted.vcf.gz

## 6. perform the final merging ##
echo -e "\nAlmost there, last step in progress..."
bcftools concat -a -Oz -o $final_out $gsa_new_prefix.sorted.vcf.gz $wes_new_prefix.vcf.gz --threads 5

if [ -f $final_out ]; then    
    echo -e "\nAll done, duration: ${SECONDS}"
else
    echo -e "\nERROR: failed to complete!"
fi
