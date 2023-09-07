#!/usr/bin/env bash
#
# @description phase genotyping array calls using SHAPEIT5
#
# GK - Jul 14th, 2023
# based on https://github.com/frhl/wes_ko_ukbb/blob/main/scripts/phasing/phasing/01_phase_common.sh
# genetic maps downloaded from https://github.com/odelaneau/shapeit5/tree/main/resources/maps
# see README.md for more details

# ISSUE: have to be careful with the suffix as vcf, vcf.gz, or '' cause Shapeit5 to crash!
#        see https://github.com/odelaneau/shapeit5/issues/48
# ISSUE: with missing AC field; should use `bcftools +fill-tags IN -Oz -o OUT -- -t AN,AC` to update - SOLVED
# ISSUE: [W::vcf_parse_filter] FILTER 'ExcessHet' is not defined in the header - SOLVED

SHAPEIT_phase_common='/software/team281/bin/shapeit5/phase_common_static'
module load common-apps/bcftools/1.16
export BCFTOOLS_PLUGINS='/software/team281/bin/bcftools_plugins'

tranche="39k"
threads=10
if [[ -z $1 ]]; then
    chr=$LSB_JOBINDEX
else
    chr=$1
fi

readonly in_file="merged_genotypes/chr$chr.merged_WES_GSA.vcf.gz"
readonly out_dir="/FIXTHIS/phasing"
readonly out_prefix="${out_dir}/phased_genotypes_common/GNH_${tranche}.chr${chr}.phased_common"
readonly log="${out_prefix}.log"
readonly gmap="/FIXTHIS/phasing/genetic_maps/chr${chr}.b38.gmap.gz"
readonly phased="${out_prefix}.bcf"

pbwt_modulo=0.1 # deault is 0.1, in shapeit4 with sequencing it is 0.0002
pbwt_depth=4 # deault is 4
pbwt_mac=5 # deafult is 5
pbwt_mdr=0.1 # default is 0.1
min_maf=0.001 # no default value

mkdir -p ${out_dir}/phased_genotypes_common

if [ ! -f ${in_file}.csi ]; then
  # WARNING: should first call 00_merge_wes_array.sh to merge WES with GSA array data)
  SECONDS=0
  echo "Creating a new VCF for chr-$chr with the suitable AC field..."
  bcftools view $in_file | bcftools +fill-tags -Oz -o $in_file.tmp -- -t AN,AC
  mv $in_file.tmp $in_file
  echo "Generating the index for $in_file..."
  bcftools index $in_file

  # echo "[deprecated] Creating a new VCF for chr-$chr with the suitable INFO field..."
  # bcftools +fill-tags ${in_dir}/chr${chr}_hard_filters.vcf.gz -r chr22 15690598 20690598  -Oz -o $in_file -- -t AN,AC
  # bcftools view $in_file | bcftools +fill-tags -Oz -o $in_file.tmp -- -t AN,AC
  # mv $in_file.tmp $in_file
  echo "Done: ${SECONDS}"
fi

if [ -f $out_dir/GNH_100trios.parents ]; then 
  tag=no100trios

  in_file_new="merged_genotypes/chr$chr.merged_WES_GSA.$tag.bcf"
  if [ ! -f $in_file_new ]; then
    echo "Generating a new BCF for $tag..."
    bcftools view $in_file -S ^$out_dir/GNH_100trios.parents -Ob -o $in_file_new
    bcftools index $in_file_new
    if [ ! -f $work_dir/samples.WES.final ]; then
      bcftools query $in_file_new --list samples > $work_dir/samples.WES.final
    fi
  fi

  SECONDS=0
  ${SHAPEIT_phase_common} \
    --input ${in_file_new} \
    --map ${gmap} \
    --region "chr${chr}" \
    --thread ${threads} \
    --output "${out_prefix}.$tag.bcf" \
    --log "${out_prefix}.$tag.log" \
    --pbwt-modulo ${pbwt_modulo} \
    --pbwt-depth ${pbwt_depth} \
    --pbwt-mac ${pbwt_mac} \
    --pbwt-mdr ${pbwt_mdr} \
    --filter-maf ${min_maf} \
    && echo "Finished phasing variants for chr${chr}, out: ${phased}, duration: ${SECONDS}" \
    || echo "ERROR: Phasing variants failed for chr${chr}, duration: ${SECONDS}"
    module purge

else

  echo "Phasing including all samples available."
  if [ ! -f ${phased} ]; then
    SECONDS=0
    # set_up_shapeit5
    ${SHAPEIT_phase_common} \
      --input ${in_file} \
      --map ${gmap} \
      --region "chr${chr}" \
      --thread ${threads} \
      --output ${phased} \
      --pbwt-modulo ${pbwt_modulo} \
      --pbwt-depth ${pbwt_depth} \
      --pbwt-mac ${pbwt_mac} \
      --pbwt-mdr ${pbwt_mdr} \
      --filter-maf ${min_maf} \
      --log "${log}" \
      && echo "Finished phasing variants for chr${chr}, out: ${phased}, duration: ${SECONDS}" \
      || echo "ERROR: Phasing variants failed for chr${chr}, duration: ${SECONDS}"
      module purge
  else
    echo "Nothing to phase, $phased already exists."
  fi

fi

# if [ ! -f "${out}.tbi" ]; then
#   module purge
#   module load BCFtools/1.12-GCC-10.3.0
#   make_tabix "${out}" "tbi"
# fi

# run as:
# bsub -J phase_common[20-22] -o logs/phase_common.%I.%J -q normal -R"select[mem>8000] rusage[mem=8000]" -M8000 -n 10 bash 01_phase_common.sh


