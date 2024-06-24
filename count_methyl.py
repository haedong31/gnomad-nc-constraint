import hail as hl
from pathlib import Path

def main():
    data_dir = Path('/mnt/isilon/immgen_res/data/gnomad-nc-constraint-v31-paper')
    hl.init()

    context_ht = hl.read_table(f'{data_dir}/context_downsampled_1000.ht')
    grouping = hl.struct(locus=context_ht.locus, context=context_ht.context, ref=context_ht.ref, 
                        alt=context_ht.alt, methyl_level=context_ht.methyl_level)
    output = {'variant_count': hl.agg.count()}

    pos_methyl_ht = context_ht.group_by(**grouping).aggregate(**output)
    pos_methyl_ht.write('data/pos_methyl_downsampled_1000.ht')

if __name__ == '__main__':
    main()
