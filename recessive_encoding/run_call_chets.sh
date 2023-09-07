#!/bin/bash

CHR=$1
work_dir="/FIXTHIS/phasing"
cpp_dir="/software/team281/gk18/call_chets/"
out_dir="$work_dir/recessive_encoding/call_chets"

BCF="$work_dir/phased_genotypes_rare/GNH_39k.chr$CHR.phased_rare.no100trios.bcf"
annot="/FIXTHIS/annotations/brava_GNH_44k_chr$CHR.txt"
genotypes="$out_dir/chr$CHR.nonref.txt.gz"

if [ ! -f $genotypes ]; then 
    echo "Generating a file with all genotypes..."
    # bcftools query $BCF --list samples > $out_dir/samples.txt
    # bcftools view $BCF --max-af 0.01 -Ou | bcftools query -i'GT="alt"' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT\n]' | gzip > $genotypes

    bcftools view --max-af 0.01 $BCF -Oz -o $out_dir/chr$CHR.phased.vcf.gz
    $cpp_dir/get_non_ref_sites $out_dir/chr$CHR.phased.vcf.gz $genotypes
fi

# grep var brava_GNH_44k_chr22.txt | awk '{print $1}' > example.genes.txt
# rm gene_map.txt
# cat example.genes.txt | while read gene; do grep $gene $annot| head -1 | tr ' ' '\n' | awk -v g=$gene 'NR>2{print $1, g}' >> gene_map.txt; done

for consq in pLoF damaging_missense synonymous; do
    python prepare_genemap.py -a $annot -c $consq -o $out_dir/gene_map.chr$CHR.$consq.txt
    $cpp_dir/call_chets -g $genotypes -m $out_dir/gene_map.chr$CHR.$consq.txt > $out_dir/chr$CHR.gen_all.$consq.txt
    for bial in cis chet hom; do
        grep $bial $out_dir/chr$CHR.gen_all.$consq.txt > $out_dir/chr$CHR.$bial.$consq.txt
        tmp=$(wc -l $out_dir/chr$CHR.$bial.$consq.txt | cut -d' ' -f1)
        echo "$bial-$consq events found: $tmp"
    done
    $cpp_dir/encode_vcf -i $out_dir/chr$CHR.gen_all.$consq.txt -s call_chets/samples.txt -m additive  | bgzip > $out_dir/chr$CHR.additive.$consq.vcf.gz
    $cpp_dir/encode_vcf -i $out_dir/chr$CHR.gen_all.$consq.txt -s call_chets/samples.txt -m recessive | bgzip > $out_dir/chr$CHR.recessive.$consq.vcf.gz
    # rm $out_dir/chr$CHR.gen_all.$consq.txt
done