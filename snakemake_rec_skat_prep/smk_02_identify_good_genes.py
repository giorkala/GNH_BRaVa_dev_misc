"""
This script will survey each gene and identify those with enough carriers to be processed in the next step.
For each of these genes, two lists will be generated: one with biallelic carriers and one with the rest.
"""
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Generate a list of genes with enough biallelic carriers to process.")
parser.add_argument("--geno", "-g", help="file with biallelic genotypes (output of `call_chets`)", required=True, type=str)
parser.add_argument("--min_carriers", "-m", help="minimum number of carriers", required=True, type=int)
parser.add_argument("--samples", "-w", help="list of available samples", required=True, type=str)
parser.add_argument("--out", "-o", help="filename for the output", required=True, type=str)
parser.add_argument("--list_carriers", help="prefix to save lists of carriers per good gene (optional)", type=str, default=None)
# parser.add_argument("--chrom", "-c", help="chromosome to work with", required=True, type=str)
args = parser.parse_args()

df = pd.read_csv( args.geno, sep='\t', header=None, names=['ID','chr','gene','gt', 'ac'])
samples_all = pd.read_csv(args.samples, sep='\t', header=None)[0].values
# restrict to avail samples (in case of sub-sampling)
df = df.query('ID in @samples_all')
df = df.query('gt=="hom" or gt=="chet"')
assert len(df)>0, "No samples left!"

# genes_all = df.gene.unique()
tmp = df.gene.value_counts()
bad_genes = tmp[ tmp < args.min_carriers].index
good_genes = tmp[ tmp >= args.min_carriers].index

print(f"Found {len(good_genes)} genes with at least {args.min_carriers} carriers.")

# make lists of good/bad genes for the current chromosome
df.query('gene in @good_genes')[['chr','gene']].drop_duplicates().to_csv(args.out, sep='\t', header=False, index=False)
df.query('gene in @bad_genes' )[['chr','gene']].drop_duplicates().to_csv(args.out+'.skip', sep='\t', header=False, index=False)
# with open(args.out, 'w') as f:
#     for gene in good_genes:
#         f.write(f"{args.chrom}\t{gene}\n")

# now make lists of biallelic/other carriers for each good gene
if args.list_carriers is not None:
    print(f"Will save lists of carriers in {args.list_carriers}/...")
    for gene in good_genes:
        samples_biallelic = df.query(f'gene=="{gene}"').ID.values
        samples_other = list( set(samples_all).difference(samples_biallelic) )

        with open(f"{args.list_carriers}/samples.{gene}.biallelic", 'w') as f:
        # with open(f"{args.list_carriers}/samples.{genes_all.loc[gene].chr}.{gene}.biallelic", 'w') as f:
            for carrier in samples_biallelic:
                f.write(f"{carrier}\n")

        with open(f"{args.list_carriers}/samples.{gene}.other", 'w') as f:
        # with open(f"{args.list_carriers}/samples.{genes_all.loc[gene].chr}.{gene}.other", 'w') as f:
            for carrier in samples_other:
                f.write(f"{carrier}\n")

# end-of-script