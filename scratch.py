import hail as hl
from pathlib import Path

data_dir = Path('./data')

mt = hl.read_table(str(data_dir/'context_downsampled_1000.ht'))
