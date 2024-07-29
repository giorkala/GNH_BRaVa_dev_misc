
# Snakemake pipeline for preparing genotypes to enable Recessive SKAT analysis in Regenie

The pipeline consists of 4 small steps, the third one being the longest one:
    1) Use `call_chets` to generate a list of biallelic individuals, once per chromosome
    2) given that list, assess which genes have enough carriers based on `min_AC`
    3) generate one VCF per gene with the prep genotypes for those genes
    4) merge all these VCFs into a single one, per chromosome, and do some retouching
Finally, concatenate all these VCFs into an exome-wide one and generate a BGEN for Regenie; this needs to be done separately. 

### Required input ####
1. BCF files with phased genotypes
2. annotation (one line per variant; can start from the SAIGE format, and convert with rule `generate_annotation_regenie_like`)
3. variant consequences to work with (e.g. pLOF + damaging_missense + other_missense)
3. Parameters: MIN_AC, max_AF, WORK_DIR

### Notes
* The rule `run_all` contains several triggering commands, one for each step. Modify accordingly to run specific parts of the pipeline.
* `generate_biallelic_carriers` needs to process the phased genotypes and extract a list of all het genotypes. This is quite slow but only needs to be performed once. Then, that needs to be passed on to `call_chets` to identify biallelic genotypes.
* `identify_good_genes` is a fairly quick step that process the output of the previous step and makes a list of genes, per chromosome, with at least four biallelic carriers to process in the next step.
* `generate_gene_vcf` is the main part, and the slowest as it needs to be invoked once for each gene. Currently, one needs to run `gene_tuples = read_genes_list(range(1,23))` before expanding the VCF filenames, meaning that we first need to run step 2, then modify `run_all`. *TODO: find a way to do this automatically*.
* `concatenate_chrom_vcf` then follows and uses bcftools to concatenate all the gene-level VCFs to one chromosome-wide. This requires a complete step3, but the current implementation doesn't check for that.
* I'm including an additional step, `generate_annotation_regenie_like`, which generates a file with the genome-wide annotation, having  one (variant,gene,consequence) tuple per line (starting from the SAIGE-like annotation file).

### how to LSF:
```
snakemake --cluster "bsub -M 4G -R 'select[mem>4G] rusage[mem=4G]' -n2 -G team281 -q normal -o run_all.stdout -e run_all.stderr" --jobs 50 --latency-wait 10 -T 2 run_all
```
this allocates two cores per task, submits up to 50 tasks simultaneously, while trying twice in case of failures.

GK - Jul 29th, 2024
"""
