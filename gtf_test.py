#Just testing gtfs
from SPORK_utils import *
gtf_path = '/scratch/PI/horence/rob/SPACHETE_dirs/SPACHETE/gtfs/hg19_gtfs'
gtfs = generate_gtfs(gtf_path)

for gtf in gtfs:
    if gtf.gene_name == 'WNT1' or gtf.gene_name == 'MLL2':
        print gtf
