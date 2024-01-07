import hail as hl
from bokeh.io import show
from pathlib import Path

hl.init()

##### Basic Table tutorial -----
hl.utils.get_movie_lens('data/')
users = hl.read_table('data/users.ht')

# Basic methods to investigate data
users.describe()
users.show()
users.count()

# Access fileds of tables with the Python attribute or index notation
users.occupation.describe()
users['occupation'].describe()
users['occupation'].show()

##### Aggregation tutorial -----
# `aggregation` provides summary quantities
# aggregators in hl.agg

users.aggregate(hl.agg.stats(users['age']))
users.aggregate(hl.agg.counter(users['occupation']))

users.filter(users['sex'] == 'M').show()
users.aggregate(hl.agg.filter(users['sex'] == 'M', hl.agg.count()))
users.aggregate(hl.agg.count_where(users['sex'] == 'M'))
users.aggregate(hl.agg.count_where((users['sex'] == 'F') & (users['occupation'] == 'scientist')))
users.aggregate(hl.agg.count_where((users['sex'] == 'F') | (users['occupation'] == 'executive')))

hist = users.aggregate(hl.agg.hist(users['age'], 10, 70, 10))
p = hl.plot.histogram(hist, legend='Age')
show(p)

##### Table joins tutorial
movies = hl.read_table('data/movies.ht')
ratings = hl.read_table('data/ratings.ht')

##### MatrixTable tutorial
hl.utils.get_1kg('data/')
mt = hl.read_matrix_table('data/1kg.mt')

##### context_downsampled_1000 -----
data_dir = Path('./data')
context_tbl = hl.read_table(str(data_dir/'context_downsampled_1000.ht'))

context_tbl.describe()
context_tbl.show()
context_tbl.count()

context_tbl.aggregate(hl.agg.collect_as_set(context_tbl['variant_type']))