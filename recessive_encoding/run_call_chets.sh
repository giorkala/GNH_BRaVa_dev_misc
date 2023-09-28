#!/bin/bash

# remember to change/select file tags

cpp_dir="/software/team281/gk18/call_chets/"
module load common-apps/bcftools/1.18

if [[ -z $1 ]]; then
    CHR=$LSB_JOBINDEX
else
    CHR=$1
fi

TAG='PP90af01'
# TAG='naive'

echo "Will create VCFs with biallelic encodings for chr-$CHR and $TAG."

work_dir="/lustre/scratch123/hgi/mdt2/projects/gnh_industry/Genes_and_Health_2023_02_44k/phasing"
out_dir="$work_dir/recessive_encoding/vcf_ready"
tmp_dir="$work_dir/recessive_encoding/sandbox"
BCF="$work_dir/phased_genotypes_rare/GNH_39k.chr$CHR.phased_rare.no100triosNew.bcf"
annot="/lustre/scratch126/humgen/teams/martin/users/gk18/work_GNH/annotations/brava_GNH_44k_chr$CHR.txt"
genotypes="$tmp_dir/chr$CHR.$TAG.txt.gz"

if [ ! -f $tmp_dir/samples.txt ]; then 
    bcftools query $BCF --list samples > $tmp_dir/samples.txt
fi

if [ ! -f $genotypes ]; then 
    echo "Generating a file with all genotypes..."
    # select all heterozygous genotypes with PP>0.90 or missing (which implies a common variant)
    bcftools view $BCF --max-af 0.01 -Ou | bcftools query -i'GT="alt"' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT %PP\n]' | awk '($4=="." || $4>0.90)' | gzip > $genotypes
    # older ways to do it:
    # bcftools view $BCF --max-af 0.01 -Ou | bcftools query -i'GT="alt"' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT\n]' | gzip > $genotypes
    # bcftools view --max-af 0.01 $BCF -Oz -o $tmp_dir/chr$CHR.phased_maxmaf01.vcf.g
    # $cpp_dir/get_non_ref_sites $tmp_dir/chr$CHR.phased_maxmaf01.vcf.gz $genotypes
fi

for consq in pLoF damaging_missense synonymous; do
    python prepare_genemap -a $annot -c $consq -o $tmp_dir/gene_map.chr$CHR.$consq.txt
    $cpp_dir/call_chets -g $genotypes -m $tmp_dir/gene_map.chr$CHR.$consq.txt > $tmp_dir/chr$CHR.gen_$TAG.$consq.txt
    $cpp_dir/encode_vcf -i $tmp_dir/chr$CHR.gen_$TAG.$consq.txt -s $tmp_dir/samples.txt -m additive  | bgzip > $out_dir/chr$CHR.additive_$TAG.$consq.vcf.gz
    $cpp_dir/encode_vcf -i $tmp_dir/chr$CHR.gen_$TAG.$consq.txt -s $tmp_dir/samples.txt -m recessive | bgzip > $out_dir/chr$CHR.recessive_$TAG.$consq.vcf.gz
done

# repeat for pLoF & damaging_missense
consq="pLoF_damaging"
python prepare_genemap -a $annot -c pLoF damaging_missense -o $tmp_dir/gene_map.chr$CHR.$consq.txt
$cpp_dir/call_chets -g $genotypes -m $tmp_dir/gene_map.chr$CHR.$consq.txt > $tmp_dir/chr$CHR.gen_$TAG.$consq.txt
$cpp_dir/encode_vcf -i $tmp_dir/chr$CHR.gen_$TAG.$consq.txt -s $tmp_dir/samples.txt -m additive  | bgzip > $out_dir/chr$CHR.additive_$TAG.$consq.vcf.gz
$cpp_dir/encode_vcf -i $tmp_dir/chr$CHR.gen_$TAG.$consq.txt -s $tmp_dir/samples.txt -m recessive | bgzip > $out_dir/chr$CHR.recessive_$TAG.$consq.vcf.gz
    
# print summaries:
for bial in chet hom cis; do
    for consq in pLoF pLoF_damaging damaging_missense synonymous; do
        tmp=$(grep $bial $tmp_dir/chr$CHR.gen_$TAG.$consq.txt | wc -l | cut -d' ' -f1)
        echo "$bial-$consq events found: $tmp"
    done
done

# rm $out_dir/chr$CHR.gen_all.*.txt

# bsub -J biall_encoding.[1-22] -o biall_encoding.%I.%J -q normal -R"select[mem>8000] rusage[mem=8000]" -M8000 -n 1 bash /nfs/users/nfs_g/gk18/SCRIPTS/GNH_BRaVa_dev_misc/recessive_encoding/run_call_chets.sh
