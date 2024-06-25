import hail as hl
import pandas as pd
from pathlib import Path

# hl.init()

def get_methyl_adj(chrom, pos, context, ref, alt, context_ht, methyl_adj_tbl=None):
    data_dir = Path('/mnt/isilon/immgen_res/data/gnomad-nc-constraint-v31-paper')
    if methyl_adj_tbl is None:
        methyl_adj_tbl = pd.read_csv(data_dir/'methyl_lookup.txt', sep='\t')
    
    locus_expr = hl.locus(chrom, pos, reference_genome='GRCh38')
    context_4input = context_ht.filter((context_ht.context == context) & (context_ht.ref == ref) & (context_ht.alt == alt) & (context_ht.locus == locus_expr))

    methyl_level = context_4input.aggregate(hl.agg.collect(context_4input.methyl_level))
    methyl_level = methyl_level[0]
    row_slice = methyl_adj_tbl[(methyl_adj_tbl['context'] == context) & (methyl_adj_tbl['ref'] == ref) & (methyl_adj_tbl['alt'] == alt) & (methyl_adj_tbl['methylation_level'] == methyl_level)]

    row_slice = row_slice[['context', 'ref', 'alt', 'methylation_level', 'methyl_factor']].reset_index()
    row_pos = pd.DataFrame({'chrom': chrom, 'pos':pos}, index=[0])
    
    return pd.concat([row_pos, row_slice], axis=1)
