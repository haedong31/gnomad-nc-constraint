import hail as hl
import pandas as pd
from pathlib import Path

def get_methyl_adj(chrom, start, end, context_ht=None, methyl_adj_tbl=None):
    data_dir = Path('/mnt/isilon/immgen_res/data/gnomad-nc-constraint-v31-paper')
    if methyl_adj_tbl is None:
        methyl_adj_tbl = pd.read_csv(data_dir/'methyl_lookup.txt', sep='\t')
        if 'methylation_level' in methyl_adj_tbl.columns:
            methyl_adj_tbl.rename(columns={'methylation_level':'methyl_level'}, inplace=True)

    if context_ht is None:
        hl.init()
        context_ht = hl.read_table(str(data_dir/'context_prefiltered.ht'))

    intervals = [hl.parse_locus_interval(f'{chrom}:{str(start)}-{str(end)}', reference_genome='GRCh38')]
    context_sub_df = hl.filter_intervals(context_ht, intervals)
    context_sub_df = context_sub_df.select(
        context_sub_df.context,
        context_sub_df.ref,
        context_sub_df.alt,
        context_sub_df.methyl_level
    ).to_pandas()
    group_vars = ['context','ref','alt','methyl_level']

    context_sub_df['locus'] = context_sub_df['locus'].astype(str)
    context_sub_df[['chrom','pos']] = context_sub_df['locus'].str.split(':', expand=True)
    context_sub_df.drop(columns='locus', inplace=True)
    context_sub_df = context_sub_df[['chrom','pos','alleles']+group_vars]
    
    return pd.merge(context_sub_df, methyl_adj_tbl[group_vars+['methyl_factor']], how='left', on=group_vars)

