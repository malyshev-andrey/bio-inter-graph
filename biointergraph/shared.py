import os
from pathlib import Path
from typing import Callable, IO

from joblib import Memory
import pandas as pd
from tqdm.auto import tqdm

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
    'lnc_RNA': 'lncRNA',
    'protein_coding': 'mRNA'
}

GOOGLE_DRIVE_URL = 'https://drive.usercontent.google.com/download?id={id}&export=download&confirm=t'

CHUNKSIZE = 10**4


# cache config
cache_dir = os.path.join(
    os.getenv('XDG_CACHE_HOME', os.path.expanduser('~/.cache')),
    'bio-inter-graph'
)
os.makedirs(cache_dir, exist_ok=True)
memory = Memory(cache_dir, verbose=0)


def _read_tsv(
        filepath_or_buffer: str | Path | IO[str], *,
        filter_func: Callable[[pd.DataFrame], pd.DataFrame] = lambda df: df,
        chunksize: int | None = CHUNKSIZE,
        desc: str | None = None,
        **kwargs
    ) -> pd.DataFrame:

    read_csv_kwargs = dict(
        sep='\t',
        dtype='str'
    )
    read_csv_kwargs.update(kwargs)

    if chunksize is None:
        return filter_func(pd.read_csv(filepath_or_buffer, **read_csv_kwargs))

    if desc is None:
        desc = filepath_or_buffer if isinstance(filepath_or_buffer, str) else ''
    with tqdm(desc=desc, unit='row') as progress_bar:
        result = []
        for chunk in pd.read_csv(filepath_or_buffer, chunksize=chunksize, **read_csv_kwargs):
            progress_bar.update(chunk.shape[0])
            result.append(filter_func(chunk))
        result = pd.concat(result)
    return result
