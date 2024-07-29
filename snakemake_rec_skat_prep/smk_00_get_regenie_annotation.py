"""
This program converts a gene-variant annotation from the SAIGE-BRaVa format to the Regenie-like format,
having one line per (variant, gene, consequence) tuple. 
It also filters out variants with non-coding consequences, or any that is not passed by the user.
"""
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Generate a list of (variant, gene, consequence) for subsequent analyses.")
parser.add_argument("--anno_prefix", "-a", help="path to chrom-based annotation", required=True, type=str)
parser.add_argument("--consq", "-c", nargs='+', help="class of consequence to process (list if many)", required=True, type=str, default=None)
parser.add_argument("--out", "-o", help="filename for the output", required=True, type=str)
parser.add_argument("--chroms", help="number of chromosomes to process", type=int, default=0)
parser.add_argument("--drop_duplicates", help="drop variants in multiple genes by keeping only the first instance.", action='store_true', default=False)

args = parser.parse_args()

if args.chroms==0:
    chroms_all = range(1,23)
    print("Will work for all 22 chromosomes.")
elif args.chroms in range(1,23):
    chroms_all = [args.chroms]
    print(f"Will only work for chromosome {args.chroms}.")
else:
    print("Invalid number of chromosomes. Exiting...")
    exit(1)

df_snp_gene = {}

for c in chroms_all:
    CHR='chr'+str(c)
    print(c) #, end=', ')
    with open( f'{args.anno_prefix}.{CHR}.saige.txt', 'r' ) as fin:

        for line in fin:
            tokens = line.split()
            if tokens[1]=='var' :
                gene_set_vars = tokens[2:]
            else:
                gene_set_anns = tokens[2:]

                df_snp_gene[tokens[0]] = pd.DataFrame( {"SNPID":gene_set_vars, "Class": gene_set_anns, "Gene":tokens[0], "chr":CHR } )
        
df_snp_annotated = pd.concat( df_snp_gene, ignore_index=True )
print('\nNumber of tuples loaded:', len(df_snp_annotated))

print("Dropping variants not in", *args.consq)
df_snp_annotated = df_snp_annotated.query('Class in @args.consq')
print("Number of variants after filtering:", len(df_snp_annotated))

if args.drop_duplicates:
    m0 = len(df_snp_annotated)
    df_snp_annotated = df_snp_annotated.drop_duplicates( subset='SNPID' )
    print(f"Removed {m0-len(df_snp_annotated)} duplicates.")
else:
    print("Final number of annotated variants:", len(df_snp_annotated))

df_snp_annotated[['SNPID','Gene','Class']].to_csv(args.out, sep=' ', index=False, header=None)
print("Output written to", args.out, 'All done!')
# end-of-script