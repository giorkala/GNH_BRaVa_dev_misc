#!/bin/bash
#
# prepape a file with genotypes for trios and phase it for SER analysis
# no need to discriminate between rare/common here, as the lowest maf will be >1/600
# input: raw_vcf, array genotypes, list of trios
# + samples.update_ids_{gsa,wes}.txt (a file to change sample IDs, optional)
# + samples.WES_in_GSA.fam (an index of WES samples who are also genotyped, optional)

module load common-apps/bcftools/1.16
plink=/software/team281/bin/plink
mem=16000

### input ###
chr=$1
tag=$2
work_dir=$3
pedigree="$work_dir/$tag.100trios.pedigree"

input_wes="FIXTHIS/filtered_vcfs/chr${chr}_hard_filters.vcf.gz"
input_gsa="FIXTHIS/GSA_QCGNH_GSAv3EAMD_DupExcl_autosome_atcgSNP"
final_out="$work_dir/sandbox/$tag.trios.prepared.chr$chr.bcf"

wes_new_prefix="$work_dir/sandbox/$tag.trios.chr$chr.exome"
gsa_new_prefix="$work_dir/sandbox/$tag.trios.chr$chr.array"

# first, generate a few lists of IDs
if [ ! -f $work_dir/$tag.100trios.parents ]; then 
    echo "Generate a list of $tag and the corresponding parents"
    # NOTE: This is now performed with python according to mendel-errors
    cut -f2,3 $pedigree | tr "\t" "\n" |  sort | uniq > $work_dir/$tag.100trios.parents 
    cat $pedigree | tr "\t" "\n" | sort  | uniq > $work_dir/$tag.100trios.samples
fi
echo "Samples to keep for trio-phasing:" $(wc -l $work_dir/$tag.100trios.samples | awk '{print $1}')

# second, make a new VCF focusing on the new set of samples and MAC>0 ##
if [ ! -f $wes_new_prefix.vcf.gz ]; then
    echo -e "\nCalling BCFtools to prepare the WES input."
    bcftools view $input_wes -S $work_dir/$tag.100trios.samples -Ou | bcftools \
    reheader --samples $work_dir/samples.update_ids_wes.txt | bcftools \
    filter --exclude 'MAF==0.0' -Ob -o $wes_new_prefix.vcf.gz
fi

## third, prepare the array genotypes
# focusing on the new set of samples (on the original index) - but extracting all samples (not just the trios)
# now subselect individuals, change from the original index to the VCF one, replace "22" with "chr22", and make the vcf in one go
$plink --bfile $input_gsa --chr $chr \
    --maf 0.001 \
    --update-ids $work_dir/samples.update_ids_gsa.txt \
    --keep $work_dir/samples.WES_in_GSA.fam \
    --real-ref-alleles \
    --recode vcf-fid bgz \
    --output-chr chr26 \
    --out $gsa_new_prefix \
    --threads 5 --memory $mem

# fourth, make a list of overlapping variants to exclude from GSA
bcftools view -HG $wes_new_prefix.vcf.gz | awk '{ print $1":"$2}' > $wes_new_prefix.snpid
bcftools view -HG $gsa_new_prefix.vcf.gz | awk '{ print $1":"$2}' > $gsa_new_prefix.snpid
comm -23 <( sort $gsa_new_prefix.snpid) <( sort $wes_new_prefix.snpid) > $gsa_new_prefix.snpid_unique
sed -i 's/:/\t/g' $gsa_new_prefix.snpid_unique
a=$(wc -l $gsa_new_prefix.snpid | awk '{print $1}')
b=$(wc -l $gsa_new_prefix.snpid_unique | awk '{print $1}')
echo "Overlapping variants to remove:" $((a-b))

# finally, merge the two VCFs 
echo -e "\nProceeding with BCFtools and final merging."
bcftools index $gsa_new_prefix.vcf.gz
bcftools index $wes_new_prefix.vcf.gz
# sort individuals in GSA's VCF according to the former
bcftools view $gsa_new_prefix.vcf.gz -S $work_dir/$tag.100trios.samples -R $gsa_new_prefix.snpid_unique -Oz -o $gsa_new_prefix.sorted.vcf.gz
bcftools index $gsa_new_prefix.sorted.vcf.gz
echo -e "\nAlmost there, last step in progress..."
bcftools concat -a -Oz -o $final_out.nosnpID $gsa_new_prefix.sorted.vcf.gz $wes_new_prefix.vcf.gz --threads 5

# update the SNP IDs to a uniform index
bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%ALT' $final_out.nosnpID -Ob -o $final_out
# mv $final_out.newSNPID.bcf $final_out
bcftools index -f $final_out
