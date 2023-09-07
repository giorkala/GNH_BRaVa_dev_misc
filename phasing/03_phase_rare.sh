#!/usr/bin/env bash
#
# @description phase genotyping array calls using SHAPEIT5
#
# GK - Jul 21st, 2023
# based on https://github.com/frhl/wes_ko_ukbb/blob/main/scripts/phasing/phasing/02_phase_chunks
#   and on https://odelaneau.github.io/shapeit5/docs/documentation/phase_rare/

# genetic maps downloaded from https://github.com/odelaneau/shapeit5/tree/main/resources/maps
#  chunk lists downloaded from https://github.com/odelaneau/shapeit5/tree/main/resources/chunks/b38/4cM
# see README.md for more details

SHAPEIT_rare='/software/team281/bin/shapeit5/phase_rare_static_2c9e551'
module load common-apps/bcftools/1.16
threads=5

tranche="39k"
tag="no100trios"
if [[ -z $1 ]]; then
    chr=$LSB_JOBINDEX
else
    chr=$1
fi

## paramters for phasing with shapeit
phased_set_error="0.0001" # this is not needed 
pbwt_min_mac=2 # for shapeit5r
pbwt_depth=4 # 5
pbwt_modulo=0.1 # default is 0.1 but 0.0004 ( 0.2 / 50 ) is default value when using --sequencing arugment
pbwt_mdr=0.1
pop_effective_size=15000

work_dir="/FIXTHIS/phasing"
gmap="${work_dir}/genetic_maps/chr${chr}.b38.gmap.gz"

# vcf to phase
vcf_dir="/FIXTHIS/filtered_vcfs"
vcf_raw="${vcf_dir}/chr${chr}_hard_filters.vcf.gz"
bcf_to_phase="$work_dir/sandbox/chr${chr}.rare_prepared.$tag.bcf" #vcf.gz"
# vcf to use as scaffold
phased_scaffold="$work_dir/phased_genotypes_common/GNH_${tranche}.chr${chr}.phased_common.$tag.bcf"

# Output paths
out_dir=$work_dir/phased_genotypes_rare/chunks
out_prefix="${out_dir}/GNH_${tranche}.chr${chr}.phase_rare.$tag"
final_out="$work_dir/phased_genotypes_rare/GNH_${tranche}.chr${chr}.phased_rare.$tag.bcf"

mkdir -p $out_dir
SECONDS=0

# first, we need to make a new VCF containing the rare variants and the same set of samples
if [ ! -f $bcf_to_phase ]; then
    echo "Preparing the input for rare variants: $bcf_to_phase."
    # this is based on the "final" set of samples generated by `01_phase_common.sh`
    bcftools reheader $vcf_raw --samples $work_dir/samples.update_ids_wes.txt | bcftools view -S $work_dir/samples.WES.final -Ou | bcftools filter --exclude 'MAF>0.001' -Ob -o $bcf_to_phase
    bcftools index $bcf_to_phase
    # bcftools view $gsa_new_prefix.vcf.gz -S $work_dir/samples.WES.final -R $gsa_new_prefix.snpid_unique -Oz -o $gsa_new_prefix.sorted.vcf.gz
    echo "Done with preprocessing, duration: ${SECONDS}."
    echo -e "\n################\n"
fi

# second, phase chunks of genotypes
# awk 'NR==1{print $0}'
cat $work_dir/chunks_b38_4cM/chunks_chr$chr.txt | while read LINE; do
	CHK=$(echo $LINE | awk '{ print $1; }')
	SRG=$(echo $LINE | awk '{ print $3; }')
	IRG=$(echo $LINE | awk '{ print $4; }')

    phased_chunk=$out_prefix.chunk$CHK.bcf
    echo "Proceeding with chunk $CHK."

    if [ ! -f ${phased_chunk} ]; then
        $SHAPEIT_rare \
            --input $bcf_to_phase \
            --input-region $IRG \
            --scaffold $phased_scaffold \
            --scaffold-region $SRG \
            --map $gmap \
            --pbwt-mac $pbwt_min_mac \
            --pbwt-depth-rare $pbwt_depth \
            --pbwt-modulo $pbwt_modulo \
            --pbwt-mdr $pbwt_mdr \
            --effective-size $pop_effective_size \
            --output $phased_chunk \
            --thread $threads
            # --phased-set-error $phased_set_error 
    else
        echo "Nothing to phase, $phased_chunk already exists."
    fi
done

echo "Done with phasing, duration: ${SECONDS}."

# third, concatenate the phased chunks
if [ -f ${phased_chunk} ]; then
    ls -1v $out_prefix.chunk*.bcf > $work_dir/sandbox/files.chr$chr
    bcftools concat -n -Ob -o $final_out -f $work_dir/sandbox/files.chr$chr
    bcftools index $final_out
    rm $work_dir/sandbox/files.chr$chr
fi

# as a last step, assess phasing
if [ -f $final_out ]; then
    echo -e "\nSubmitting a job for assessing phasing accuracy:"
    bsub -J assess_phasing.$chr -o logs/assess_phasing.$chr.%J -q normal -R"select[mem>4000] rusage[mem=4000]" -M4000 bash 04_assess_phasing.sh $chr
fi

# run as 
# bsub -J phase_rare[20-22] -o logs/phase_rare.%I.%J -q normal -R"select[mem>12000] rusage[mem=12000]" -M12000 -n 5 bash 03_phase_rare.sh
