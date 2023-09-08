#!/bin/bash

cpp_dir="/software/team281/gk18/call_chets/"
module load common-apps/bcftools/1.18

if [[ -z $1 ]]; then
    CHR=$LSB_JOBINDEX
else
    CHR=$1
fi
echo "Will create VCFs with biallelic encodings for chr-$chr."

work_dir="/FIXTHIS/phasing"
out_dir="$work_dir/recessive_encoding/vcf_ready"
tmp_dir="$work_dir/recessive_encoding/sandbox"
BCF="$work_dir/phased_genotypes_rare/GNH_39k.chr$CHR.phased_rare.no100trios.bcf"
annot="/FIXTHIS/annotations/brava_GNH_44k_chr$CHR.txt"
genotypes="$tmp_dir/chr$CHR.nonref.txt.gz"

if [ ! -f $genotypes ]; then 
    echo "Generating a file with all genotypes..."
    bcftools query $BCF --list samples > $tmp_dir/samples.txt
    # bcftools view $BCF --max-af 0.01 -Ou | bcftools query -i'GT="alt"' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT\n]' | gzip > $genotypes

    bcftools view --max-af 0.01 $BCF -Oz -o $tmp_dir/chr$CHR.phased_maxmaf01.vcf.gz
    $cpp_dir/get_non_ref_sites $tmp_dir/chr$CHR.phased_maxmaf01.vcf.gz $genotypes
fi

# grep var brava_GNH_44k_chr22.txt | awk '{print $1}' > example.genes.txt
# rm gene_map.txt
# cat example.genes.txt | while read gene; do grep $gene $annot| head -1 | tr ' ' '\n' | awk -v g=$gene 'NR>2{print $1, g}' >> gene_map.txt; done

for consq in pLoF damaging_missense synonymous; do
    python prepare_genemap -a $annot -c $consq -o $tmp_dir/gene_map.chr$CHR.$consq.txt
    $cpp_dir/call_chets -g $genotypes -m $tmp_dir/gene_map.chr$CHR.$consq.txt > $tmp_dir/chr$CHR.gen_all.$consq.txt
    for bial in cis chet hom; do
        grep $bial $tmp_dir/chr$CHR.gen_all.$consq.txt > $tmp_dir/chr$CHR.$bial.$consq.txt
        tmp=$(wc -l $tmp_dir/chr$CHR.$bial.$consq.txt | cut -d' ' -f1)
        echo "$bial-$consq events found: $tmp"
    done
    $cpp_dir/encode_vcf -i $tmp_dir/chr$CHR.gen_all.$consq.txt -s $tmp_dir/samples.txt -m additive  | bgzip > $out_dir/chr$CHR.additive.$consq.vcf.gz
    $cpp_dir/encode_vcf -i $tmp_dir/chr$CHR.gen_all.$consq.txt -s $tmp_dir/samples.txt -m recessive | bgzip > $out_dir/chr$CHR.recessive.$consq.vcf.gz
    # rm $out_dir/chr$CHR.gen_all.$consq.txt
done

# repeat for pLoF & damaging_missense
consq="pLoF_damaging"
python prepare_genemap -a $annot -c pLoF damaging_missense -o $tmp_dir/gene_map.chr$CHR.$consq.txt
$cpp_dir/call_chets -g $genotypes -m $tmp_dir/gene_map.chr$CHR.$consq.txt > $tmp_dir/chr$CHR.gen_all.$consq.txt
$cpp_dir/encode_vcf -i $tmp_dir/chr$CHR.gen_all.$consq.txt -s $tmp_dir/samples.txt -m additive  | bgzip > $out_dir/chr$CHR.additive.$consq.vcf.gz
$cpp_dir/encode_vcf -i $tmp_dir/chr$CHR.gen_all.$consq.txt -s $tmp_dir/samples.txt -m recessive | bgzip > $out_dir/chr$CHR.recessive.$consq.vcf.gz

# bsub -J biall_encoding.[1-22] -o biall_encoding.%I.%J -q normal -R"select[mem>4000] rusage[mem=4000]" -M4000 -n 1 bash run_call_chets.sh
