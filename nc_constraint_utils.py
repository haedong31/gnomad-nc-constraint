#!/usr/bin/env python
# coding: utf-8

from generic import *
from constraint_basics import *


def filter_to_autosomes_par(ht: Union[hl.Table, hl.MatrixTable]) -> Union[hl.Table, hl.MatrixTable]:
    return ht.filter(ht.locus.in_autosome_or_par())        


def remove_coverage_outliers(t: Union[hl.MatrixTable, hl.Table]) -> Union[hl.MatrixTable, hl.Table]:
    '''
    Keep only loci where genome coverage was between 15 and 60
    '''
    criteria = (t.coverage_mean >= 15) & (t.coverage_mean <= 60)
    return t.filter_rows(criteria) if isinstance(t, hl.MatrixTable) else t.filter(criteria)


def count_variants(ht: hl.Table,
                   count_singletons: bool = False, count_downsamplings: Optional[List[str]] = (),
                   additional_grouping: Optional[List[str]] = (), partition_hint: int = 100,
                   omit_methylation: bool = False, return_type_only: bool = False,
                   force_grouping: bool = False, singleton_expression: hl.expr.BooleanExpression = None,
                   impose_high_af_cutoff_here: bool = False) -> Union[hl.Table, Any]:
    """
    Count variants by context, ref, alt, methylation_level
    """

    grouping = hl.struct(context=ht.context, ref=ht.ref, alt=ht.alt)
    if not omit_methylation:
        grouping = grouping.annotate(methylation_level=ht.methylation_level)
    for group in additional_grouping:
        grouping = grouping.annotate(**{group: ht[group]})

    if count_singletons:
        # singleton = hl.any(lambda f: (f.meta.size() == 1) & (f.meta.get('group') == 'adj') & (f.AC[1] == 1), ht.freq)
        if singleton_expression is None:
            singleton_expression = ht.freq[0].AC == 1

    if count_downsamplings or force_grouping:
        # Slower, but more flexible (allows for downsampling agg's)
        output = {'variant_count': hl.agg.count_where(ht.freq[0].AF <= 0.001) if impose_high_af_cutoff_here else hl.agg.count()}
        for pop in count_downsamplings:
            output[f'downsampling_counts_{pop}'] = downsampling_counts_expr(ht, pop, impose_high_af_cutoff=impose_high_af_cutoff_here)
        if count_singletons:
            output['singleton_count'] = hl.agg.count_where(singleton_expression)
            for pop in count_downsamplings:
                output[f'singleton_downsampling_counts_{pop}'] = downsampling_counts_expr(ht, pop, singleton=True)
        return ht.group_by(**grouping).partition_hint(partition_hint).aggregate(**output)
    else:
        agg = {'variant_count': hl.agg.counter(grouping)}
        if count_singletons:
            agg['singleton_count'] = hl.agg.counter(hl.agg.filter(singleton_expression, grouping))

        if return_type_only:
            return agg['variant_count'].dtype
        else:
            return ht.aggregate(hl.struct(**agg))


def annotate_variant_types(t: Union[hl.MatrixTable, hl.Table],
                           heptamers: bool = False) -> Union[hl.MatrixTable, hl.Table]:
    """
    Adds cpg, transition, and variant_type, variant_type_model columns
    """
    mid_index = 3 if heptamers else 1
    transition_expr = (((t.ref == "A") & (t.alt == "G")) | ((t.ref == "G") & (t.alt == "A")) |
                       ((t.ref == "T") & (t.alt == "C")) | ((t.ref == "C") & (t.alt == "T")))
    cpg_expr = (((t.ref == "G") & (t.alt == "A") & (t.context[mid_index - 1:mid_index] == 'C')) |
                ((t.ref == "C") & (t.alt == "T") & (t.context[mid_index + 1:mid_index + 2] == 'G')))
    if isinstance(t, hl.MatrixTable):
        t = t.annotate_rows(transition=transition_expr, cpg=cpg_expr)
    else:
        t = t.annotate(transition=transition_expr, cpg=cpg_expr)
    variant_type_expr = (hl.case()
                         .when(t.cpg, 'CpG')
                         .when(t.transition, 'non-CpG transition')
                         .default('transversion'))
    variant_type_model_expr = hl.cond(t.cpg, t.context, "non-CpG")
    if isinstance(t, hl.MatrixTable):
        return t.annotate_rows(variant_type=variant_type_expr, variant_type_model=variant_type_model_expr)
    else:
        return t.annotate(variant_type=variant_type_expr, variant_type_model=variant_type_model_expr)

 
def filter_black_regions(ht: Union[hl.MatrixTable, hl.Table]) -> Union[hl.MatrixTable, hl.Table]:
    bed = hl.import_bed('gs://gnomad-nc-constraint-v31-paper/misc/blacklist_gap2.bed', reference_genome='GRCh38', skip_invalid_intervals=True)
    return ht.filter(hl.is_defined(bed[ht.locus]), keep=False)


def annotate_genome_element(ht: Union[hl.MatrixTable, hl.Table], 
							bed_path: str) -> Union[hl.MatrixTable, hl.Table]:
    bed = hl.import_bed(bed_path, reference_genome='GRCh38', skip_invalid_intervals=True)
    ht = ht.filter(hl.is_defined(bed[ht.locus]))
    return ht.annotate(element_id = bed.index(ht.locus,all_matches=True).target).explode('element_id')









#############
import numpy as np
def sem(x,n):
    # import numpy as np
    p=float(x)/n
    # return np.sqrt(p*(1-p)/float(n))
    return math.sqrt(p*(1-p)/float(n))

def read_gs_to_lines(gs_bucket, gs_file_path):
    from google.cloud import storage
    import io
    from google.cloud import storage
    client = storage.Client()
    bucket = client.get_bucket(gs_bucket)
    blob = bucket.get_blob(gs_file_path)
    return blob.download_as_string().decode("utf-8").strip().split("\n")


def read_gs_to_df(gs_bucket, gs_file_path):
    from google.cloud import storage
    import io
    import pandas as pd
    client = storage.Client()
    bucket = client.get_bucket(gs_bucket)
    blob = bucket.get_blob(gs_file_path)
    lines = blob.download_as_string().decode("utf-8").strip()
    return pd.read_csv(io.StringIO(lines), sep="\t")

















