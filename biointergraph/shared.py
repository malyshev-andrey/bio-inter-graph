import os
import hashlib
from pathlib import Path
from typing import Callable, IO

from joblib import Memory
import requests
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


class HttpFileReader:
    def __init__(self, url: str, chunk_size: int = 1024*1024, desc: str = ''):
        def reader():
            response = requests.get(url, stream=True)
            response.raise_for_status()
            size = int(response.headers['Content-Length'])
            with tqdm(desc=desc, total=size, unit='B', unit_scale=True) as progress_bar:
                for chunk in response.iter_content(chunk_size=chunk_size):
                    progress_bar.update(len(chunk))
                    yield chunk

        self.reader = reader()
        self.buffer = b''


    def read(self, n: int = -1):
        if n < 0:
            result = self.buffer + b''.join(self.reader)
            self.buffer = b''
        else:
            buffer_size = len(self.buffer)
            chunks = []
            while buffer_size < n:
                try:
                    chunk = next(self.reader)
                    chunks.append(chunk)
                    buffer_size += len(chunk)
                except StopIteration:
                    break
            self.buffer += b''.join(chunks)
            assert len(self.buffer) == buffer_size
            result = self.buffer[:n]
            self.buffer = self.buffer[n:]
        return result


def _read_tsv(
        filepath_or_buffer: str | Path | IO[str], *,
        filter_func: Callable[[pd.DataFrame], pd.DataFrame] = lambda df: df,
        chunksize: int | None = CHUNKSIZE,
        desc: str = '',
        **kwargs
    ) -> pd.DataFrame:

    read_csv_kwargs = dict(
        sep='\t',
        dtype='str'
    )
    read_csv_kwargs.update(kwargs)

    if isinstance(filepath_or_buffer, str):
        desc = desc or filepath_or_buffer.split('://')[1]
        if chunksize is not None and filepath_or_buffer.startswith('http'):
            if 'compression' not in read_csv_kwargs:
                if filepath_or_buffer.endswith('.gz'):
                    read_csv_kwargs['compression'] = 'gzip'
                elif filepath_or_buffer.endswith('.zip'):
                    read_csv_kwargs['compression'] = 'zip'
            reader = pd.read_csv(
                HttpFileReader(filepath_or_buffer, desc=desc),
                chunksize=chunksize,
                **read_csv_kwargs
            )
            reader = map(filter_func, reader)
            return pd.concat(reader)


    if chunksize is None:
        return filter_func(pd.read_csv(filepath_or_buffer, **read_csv_kwargs))

    with tqdm(desc=desc, unit='row') as progress_bar:
        result = []
        for chunk in pd.read_csv(filepath_or_buffer, chunksize=chunksize, **read_csv_kwargs):
            progress_bar.update(chunk.shape[0])
            result.append(filter_func(chunk))
        result = pd.concat(result)
    return result


def _df_hash(df: pd.DataFrame) -> str:
    result = pd.util.hash_pandas_object(df).values
    result = hashlib.sha1(result).hexdigest()
    return result
