import os
import json
import re
from io import BytesIO
import hashlib
from pathlib import Path
from typing import Callable, IO
from time import time
from urllib.parse import urlencode, urlsplit, urlunsplit, parse_qsl

from joblib import Memory
import pandas as pd
from tqdm.auto import tqdm

import fsspec
from fsspec.callbacks import TqdmCallback

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

REMOTE_PROTOCOLS = ['https://', 'http://', 'ftp://']
REMOTE_REGEX = f"^({'|'.join(REMOTE_PROTOCOLS)})"

# cache config
cache_dir = os.path.join(
    os.getenv('XDG_CACHE_HOME', os.path.expanduser('~/.cache')),
    'bio-inter-graph'
)
fsspec_cache_dir = os.path.join(cache_dir, 'fsspec')
memory = Memory(cache_dir, verbose=0)


def _shorten_url(url: str, max_len: int = 70, ellipsis: str = "...") -> str:
    if len(url) <= max_len:
        return url

    remain = max_len - len(ellipsis)

    head_len = (remain + 1) // 2
    tail_len = remain // 2

    return url[:head_len] + ellipsis + url[-tail_len:]


def _df_hash(df: pd.DataFrame) -> str:
    result = pd.util.hash_pandas_object(df).values
    result = hashlib.sha1(result).hexdigest()
    return result


def _canonicalize_url(url: str) -> str:
    parts = urlsplit(url)

    pairs = parse_qsl(
        parts.query,
        keep_blank_values=False,
        strict_parsing=True,
        encoding="utf-8",
        errors="strict",
    )

    return urlunsplit((parts.scheme, parts.netloc, parts.path, urlencode(sorted(pairs), doseq=True), parts.fragment))


def remote_file2local(
        url: str, *,
        cache_dir: str = fsspec_cache_dir,
        progress_bar: bool = True,
        **remote_opts
    ) -> str:
    parts = url.split('::')
    n_remote_parts = sum(1 for part in parts if re.match(REMOTE_REGEX, part))
    if n_remote_parts == 0:
        return url

    assert n_remote_parts == 1
    assert re.match(REMOTE_REGEX, parts[-1])

    remote_fs, remote_path = fsspec.core.url_to_fs(_canonicalize_url(parts[-1]), **remote_opts)

    sc = fsspec.filesystem(
        "simplecache",
        fs=remote_fs,
        cache_storage=cache_dir,
    )

    cache_dir = sc.storage[-1]
    local_path = os.path.join(cache_dir, sc._mapper(remote_path))

    if not os.path.exists(local_path):
        os.makedirs(cache_dir, exist_ok=True)
        cb = TqdmCallback(
            tqdm_cls=tqdm,
            tqdm_kwargs=dict(
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
                desc='LOADING: ' + _shorten_url(remote_path)
            )
        )
        start = time()
        kwargs = dict(callback=cb) if progress_bar else {}
        remote_fs.get_file(remote_path, local_path, **kwargs)
        download_time = time() - start

        metadata_path = local_path + '.meta.json'
        with open(metadata_path, 'wt') as metadata_file:
            json.dump(
                {
                    'url': parts[-1],
                    'canonical_url': remote_path,
                    'ts': start + download_time,
                    'download_time': download_time,
                    'local_path': local_path
                },
                metadata_file,
                indent=2,
                sort_keys=True
            )

    if remote_path.endswith('.gz'):
        prefix = 'gzip::'
    elif remote_path.endswith('.zip'):
        prefix = 'zip::'
    else:
        prefix = ''

    parts[-1] = prefix + f'file://{local_path}'

    new_url = '::'.join(parts)
    return new_url


def _read_tsv(
        filepath_or_buffer: str | Path | IO[str], *,
        filter_func: Callable[[pd.DataFrame], pd.DataFrame] = lambda df: df,
        chunksize: int | None = CHUNKSIZE,
        desc: str = '',
        use_cache: bool = False,
        **kwargs
    ) -> pd.DataFrame:

    read_csv_kwargs = dict(
        sep='\t',
        dtype='str'
    )
    read_csv_kwargs.update(kwargs)

    if not desc and isinstance(filepath_or_buffer, str):
        desc = 'READING: ' + _shorten_url(filepath_or_buffer)

    if use_cache and isinstance(filepath_or_buffer, str):
        filepath_or_buffer = remote_file2local(filepath_or_buffer, progress_bar=chunksize is not None)

    if chunksize is None:
        return filter_func(pd.read_csv(filepath_or_buffer, **read_csv_kwargs))

    reader = pd.read_csv(filepath_or_buffer, chunksize=chunksize, **read_csv_kwargs)

    with tqdm(desc=desc, unit='row') as progress_bar:
        result = []
        for chunk in reader:
            progress_bar.update(chunk.shape[0])
            result.append(filter_func(chunk))
        result = pd.concat(result)
    return result
