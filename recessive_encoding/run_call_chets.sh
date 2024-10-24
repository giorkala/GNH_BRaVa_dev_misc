#!/bin/bash

# remember to change/select file tags

# module load common-apps/bcftools/1.18
module load HGI/softpack/users/jw38/primedel/1
# module load HGI/softpack/users/bh18/bcftools/1.0
cpp_dir="/software/team281/gk18/call_chets/"

if [[ -z $1 ]]; then
    CHR=$LSB_JOBINDEX
else
    CHR=$1
fi

AFthres="05"
TAG="PP90af05.s50"
# TAG="BRAVA"
# TAG="gnh_flagship"

echo "Will create VCFs with biallelic encodings for chr-$CHR and $TAG."

work_dir="/lustre/scratch123/hgi/mdt2/projects/gnh_industry/Genes_and_Health_2023_02_44k"
out_dir="$work_dir/phasing/recessive_encoding/vcf_ready"
tmp_dir="$work_dir/phasing/recessive_encoding/sandbox"
BCF="$work_dir/phasing/phased_genotypes/GNH_39k_QC.notrios.phased_all.chr$CHR.bcf"
# annot="/lustre/scratch126/humgen/teams/martin/users/gk18/work_GNH/annotations/brava_GNH_annot_v2.chr$CHR.txt"
# annot="/lustre/scratch126/humgen/teams/martin/users/gk18/work_GNH/annotations/brava_GNH_annot_RCV.chr$CHR.txt"
annot="$work_dir/brava_vep_annotation/vep105_loftee/out/gnh_39k.brava_s50_nopopmax.chr$CHR.saige.txt"
# annot="${work_dir}/phasing/recessive_encoding/annotation_regenie_analysis/gnh_39k.flagship_annotation.chr$CHR.saige.txt"
genotypes="$tmp_dir/chr$CHR.PP90af05.txt.gz"

if [ ! -f $tmp_dir/samples.txt ]; then 
    bcftools query $BCF --list samples > $tmp_dir/samples.txt
fi

if [ ! -f $genotypes ]; then 
    echo "Generating a file with all genotypes..."
    # select all heterozygous genotypes with PP>0.90 or missing (which implies a common variant)
#    bcftools view $BCF --min-af 0.$AFthres -Ou | bcftools query -i'GT="alt"' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT %PP\n]' | awk '($4=="." || $4>0.90)' | gzip > $genotypes
    bcftools view $BCF --max-af 0.$AFthres -Ou | bcftools query -i'GT="alt"' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT %PP\n]' | awk '(length($2)<50)' | awk '($4=="." || $4>0.90)' | gzip > $genotypes
    # older ways to do it:
    # $cpp_dir/get_non_ref_sites $tmp_dir/chr$CHR.phased_maxmaf01.vcf.gz $genotypes
fi

for consq in pLoF damaging_missense_or_protein_altering synonymous; do
# for consq in pLOF deleterious_missense other_missense synonymous; do
    python prepare_genemap -a $annot -c $consq -o $tmp_dir/gene_map.chr$CHR.$consq.txt
    $cpp_dir/call_chets -g $genotypes -m $tmp_dir/gene_map.chr$CHR.$consq.txt > $tmp_dir/chr$CHR.gen_$TAG.$consq.txt
    # $cpp_dir/encode_vcf -i $tmp_dir/chr$CHR.gen_$TAG.$consq.txt -s $tmp_dir/samples.txt -m additive  | bgzip > $out_dir/chr$CHR.additive_$TAG.$consq.vcf.gz
    # $cpp_dir/encode_vcf -i $tmp_dir/chr$CHR.gen_$TAG.$consq.txt -s $tmp_dir/samples.txt -m recessive | bgzip > $out_dir/chr$CHR.recessive_$TAG.$consq.vcf.gz
done

# repeat for pLoF & damaging_missense_or_protein_altering
consq="pLoF_damaging"
python prepare_genemap -a $annot -c pLoF damaging_missense_or_protein_altering -o $tmp_dir/gene_map.chr$CHR.$TAG.$consq.txt
$cpp_dir/call_chets -g $genotypes -m $tmp_dir/gene_map.chr$CHR.$TAG.$consq.txt > $tmp_dir/chr$CHR.gen_$TAG.$consq.txt

consq="nonsynonymous"
python prepare_genemap -a $annot -c pLoF damaging_missense_or_protein_altering other_missense_or_protein_altering -o $tmp_dir/gene_map.chr$CHR.$TAG.$consq.txt
$cpp_dir/call_chets -g $genotypes -m $tmp_dir/gene_map.chr$CHR.$TAG.$consq.txt > $tmp_dir/chr$CHR.gen_$TAG.$consq.txt

for consq in pLoF damaging_missense_or_protein_altering pLoF_damaging nonsynonymous synonymous; do
    $cpp_dir/encode_vcf -i $tmp_dir/chr$CHR.gen_$TAG.$consq.txt -s $tmp_dir/samples.txt -m additive  | bgzip > $out_dir/chr$CHR.add_$TAG.$consq.vcf.gz
    $cpp_dir/encode_vcf -i $tmp_dir/chr$CHR.gen_$TAG.$consq.txt -s $tmp_dir/samples.txt -m dominance | bgzip > $out_dir/chr$CHR.dom_$TAG.$consq.vcf.gz
    $cpp_dir/encode_vcf -i $tmp_dir/chr$CHR.gen_$TAG.$consq.txt -s $tmp_dir/samples.txt -m recessive | bgzip > $out_dir/chr$CHR.rec_$TAG.$consq.vcf.gz
done

# # print summaries:
# for bial in chet hom cis; do
#     for consq in pLoF pLoF_damaging damaging_missense_or_protein_altering synonymous; do
#         tmp=$(grep $bial $tmp_dir/chr$CHR.gen_$TAG.$consq.txt | wc -l | cut -d' ' -f1)
#         echo "$bial-$consq events found: $tmp"
#     done
# done

# rm $out_dir/chr$CHR.gen_all.*.txt

# bsub -J biall_encoding.[1-22] -o biallenc.%I.%J -q normal -R"select[mem>1000] rusage[mem=1000]" -M1000 -n 1 bash /nfs/users/nfs_g/gk18/SCRIPTS/GNH_BRaVa_dev_misc/recessive_encoding/run_call_chets.sh
