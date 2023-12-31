#!/bin/bash

# GK - Jul 28th, 2023
#
# based on https://odelaneau.github.io/shapeit5/docs/documentation/phase_common/

# In brief, this generates a list of trios (and sample lists for parents etc), 
# then creates BCF files for that set and phases, both for common and rare variants.
# see README.md for more details

if [[ -z $1 ]]; then
    chr=$LSB_JOBINDEX
else
    chr=$1
fi

module load common-apps/bcftools/1.16
SHAPEIT_phase_common='/software/team281/bin/shapeit5/phase_common_static_2c9e551'
SHAPEIT_rare='/software/team281/bin/shapeit5/phase_rare_static_2c9e551'
threads=5

work_dir="/FIXTHIS/phasing"
gmap="${work_dir}/genetic_maps/chr${chr}.b38.gmap.gz"
tag="100trios" # remember to change the awk command at line 34
tranche="39k"

echo -e "########################\n  Phasing for $tag \n  Working for chrom-$chr \n########################\n"

if [ ! -f $work_dir/GNH_$tag.pedigree ]; then 
    echo "Generate a list of $tag and the corresponding parents"
    # awk 'NR<101{print $2,$3,$4}' GH_44k_668-trios_QCed.mercury.consistent.fam | sed -e 's/\_2//g' -e 's/\_3//g' -e 's/\_4//g' > GNH_$tag.pedigree
    join -j 1 <(cut -f2-4 $work_dir/GNH.668trios.mercury.consistent.fam | sort) <( awk '{print $1}' $work_dir/samples.WES_in_GSA.fam | sort ) > $work_dir/temp.1
    join -1 2 -2 1 <(cut -f2-4 $work_dir/temp.1 | sort -k2) <( awk '{print $1}' $work_dir/samples.WES_in_GSA.fam | sort ) > $work_dir/temp.2
    join -1 3 -2 1 <(cut -f2-4 $work_dir/temp.2 | sort -k3) <( awk '{print $1}' $work_dir/samples.WES_in_GSA.fam | sort ) | awk 'NR<101{print $3,$2,$1}' OFS="\t" > $work_dir/GNH_$tag.pedigree 
    cut -f2,3 $work_dir/GNH_$tag.pedigree | tr "\t" "\n" |  sort | uniq > $work_dir/GNH_$tag.parents 
    cat $work_dir/GNH_$tag.pedigree | tr "\t" "\n" | sort  | uniq > $work_dir/GNH_$tag.samples
    rm $work_dir/temp.1 $work_dir/temp.2
fi

echo "Samples to keep for analysis:" $(wc -l GNH_100trios.samples | awk '{print $1}')

echo -e "\n################\n"
## phase common variants ## 


in_merged="${work_dir}/merged_genotypes/chr$chr.merged_WES_GSA.vcf.gz"
to_phase="${work_dir}/merged_genotypes/chr$chr.merged_WES_GSA.$tag.bcf"
out_prefix="${work_dir}/phased_genotypes_common/GNH_$tranche.chr${chr}.phased_common.$tag.bcf"

if [ ! -f $to_phase ]; then
    echo "Generating a new BCF for trios samples..."
    SECONDS=0
    # bcftools view $in_file -S ^GNH_$tag.parents -Ob -o chr$chr.merged.notrios.bcf
    bcftools view $in_merged -S $work_dir/GNH_$tag.samples -Ob -o $to_phase
    bcftools index $to_phase
    echo "Done with preprocessing common, duration: ${SECONDS}."
else
    echo "Input is already prepared: $to_phase."
fi

if [ ! -f $out_prefix ]; then

    pbwt_modulo=0.1 # deault is 0.1, in shapeit4 with sequencing it is 0.0002
    pbwt_depth=4 # deault is 4
    pbwt_mac=5 # deafult is 5
    pbwt_mdr=0.1 # default is 0.1
    min_maf=0.001 # no default value

    ${SHAPEIT_phase_common} \
        --input $to_phase \
        --map ${gmap} \
        --region "chr${chr}" \
        --pedigree $work_dir/GNH_$tag.pedigree \
        --thread ${threads} \
        --output ${out_prefix} \
        --log $out_prefix.log \
        --pbwt-modulo ${pbwt_modulo} \
        --pbwt-depth ${pbwt_depth} \
        --pbwt-mac ${pbwt_mac} \
        --pbwt-mdr ${pbwt_mdr} \
        --filter-maf ${min_maf}
else
    echo "Nothing to do for common, $out_prefix already exists!"
fi

echo -e "\n################\n"
## phase rare variants ##

vcf_dir="/FIXTHIS/filtered_vcfs"
vcf_raw="${vcf_dir}/chr${chr}_hard_filters.vcf.gz"
scaffold="$work_dir/phased_genotypes_common/GNH_${tranche}.chr${chr}.phased_common.$tag.bcf"
to_phase="$work_dir/sandbox/GNH_${tranche}.chr${chr}.phased_rare.$tag.bcf"

out_dir=$work_dir/phased_genotypes_rare/chunks
out_prefix="${out_dir}/GNH_${tranche}.chr${chr}.phase_rare.$tag"
final_out="$work_dir/phased_genotypes_rare/GNH_${tranche}.chr${chr}.phased_rare.$tag.bcf"

mkdir -p $out_dir

if [ ! -f $to_phase ]; then
    echo "Preparing the input for rare variants: $to_phase."
    SECONDS=0
    bcftools reheader $vcf_raw --samples $work_dir/samples.update_ids_wes.txt | bcftools view -S $work_dir/GNH_$tag.samples -Ou | bcftools filter --exclude 'MAF>0.001' -Ob -o $to_phase
    bcftools index $to_phase
    echo "Done with preprocessing rare, duration: ${SECONDS}."
else
    echo "Input is already prepared: $to_phase."
fi

if [ ! -f $final_out ]; then
    SECONDS=0

    # params for rare variants - following Lassen et al.
    pbwt_min_mac=2 # for shapeit5r
    pbwt_depth=4 # 5
    pbwt_modulo=0.1 # default is 0.1 but 0.0004 ( 0.2 / 50 ) is default value when using --sequencing arugment
    pbwt_mdr=0.1
    pop_effective_size=15000

    cat $work_dir/chunks_b38_4cM/chunks_chr$chr.txt | while read LINE; do
        CHK=$(echo $LINE | awk '{ print $1; }')
        SRG=$(echo $LINE | awk '{ print $3; }')
        IRG=$(echo $LINE | awk '{ print $4; }')

        phased_chunk=$out_prefix.chunk$CHK.bcf
        echo "Proceeding with chunk $CHK."

        if [ ! -f ${phased_chunk} ]; then
            $SHAPEIT_rare \
                --input $to_phase \
                --input-region $IRG \
                --scaffold $scaffold \
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
                # --no-index \
                # --log $out_prefix.log \
        else
            echo "Nothing to phase, $phased_chunk already exists."
        fi
    done
    echo "Done with phasing, duration: ${SECONDS}."

    # finally, concatenate the phased chunks
    if [ -f ${phased_chunk} ]; then
        ls -1v $out_prefix.chunk*.bcf > $work_dir/sandbox/files.chr$chr
        bcftools concat -n -Ob -o $final_out -f $work_dir/sandbox/files.chr$chr
        bcftools index $final_out
        rm $work_dir/sandbox/files.chr$chr
    fi

else
    echo "Nothing to do for rare, $final_out already exists!"
fi

# run as:
# bsub -J phase_related[20-22] -o logs/phase_related.%I.%J -q normal -R"select[mem>8000] rusage[mem=8000]" -M8000 -n 5 bash 010_phase_related.sh
