import os
from joblib import Memory

# schemas
GFF_COLUMNS = [
    'chr', 'source', 'type',
    'start', 'end', 'score',
    'strand', 'phase', 'attributes'
]
BED_COLUMNS = [
    'chr', 'start', 'end',
    'name', 'score', 'strand'
]
ID_TYPES = [
    'entrezgene_id',
    'ensembl_gene_id',
    'ensembl_transcript_id',
    'refseq_transcript_id'
]

UNIFY_BIOTYPES = {
    'protein-coding': 'mRNA',
    'pseudo': 'pseudogene',
    'lnc_RNA': 'lncRNA'
}

GOOGLE_DRIVE_URL = 'https://drive.usercontent.google.com/download?id={id}&export=download&confirm=t'

CHUNKSIZE = 10000


# cache config
cache_dir = os.path.join(
    os.getenv('XDG_CACHE_HOME', os.path.expanduser('~/.cache')),
    'bio-inter-graph'
)
os.makedirs(cache_dir, exist_ok=True)
memory = Memory(cache_dir, verbose=0)
