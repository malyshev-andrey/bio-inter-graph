import pandas as pd

from ..shared import memory, CHUNKSIZE
from .main import summarize_pairwise
from ..ids_mapping import id2yagid
from .karr_seq_shared import _retrieve_karr_seq_metadata, _load_single_karr_seq


def _load_karr_seq_data(
        cell_line: str|None = None, *,
        chunksize: int|None = CHUNKSIZE,
        verbose: bool = True
    ) -> pd.DataFrame:

    metadata = _retrieve_karr_seq_metadata(cell_line=cell_line)

    data = []
    for _, row in metadata.iterrows():
        if verbose:
            print(row.cell_line, row.dendrimers, row.frac, row.repl, sep='-')
            print(f'\tDownloading KARR-seq data from: {row.url}')

        data.append(_load_single_karr_seq(
            row.url,
            filter_func=lambda df: df[df['seqid1'] != df['seqid2']],
            chunksize=chunksize
        ))
        data[-1]['dendrimers'] = row.dendrimers
        data[-1]['repl'] = row.repl
        data[-1]['frac'] = row.frac
        data[-1]['cell_line'] = row.cell_line

        if verbose:
            print(f'\tShape: {data[-1].shape}')

    data = pd.concat(data)

    if verbose:
        print(f'Shape (total): {data.shape}')

    return data


def _karr_seq_data2pairwise(data: pd.DataFrame) -> pd.DataFrame:
    swap = ['seqid', 'pos', 'strand']
    swap = {f'{c}1': f'{c}2' for c in swap} | {f'{c}2': f'{c}1' for c in swap}
    data = pd.concat([
        data,
        data.rename(columns=swap)
    ])
    data = data[data['seqid1'] < data['seqid2']]
    data = data.groupby(
        ['seqid1', 'seqid2', 'repl', 'frac', 'cell_line'],
        as_index=False
    ).agg(
        n=('readID', 'size')
    )
    return data


@memory.cache
def load_karr_seq_data(pvalue: float|None):
    data = _load_karr_seq_data(cell_line='K562')
    data = _karr_seq_data2pairwise(data)
    data = data.groupby(
        ['repl', 'frac'],
        as_index=False
    ).apply(
        lambda df: summarize_pairwise(df[['seqid1', 'seqid2']], symmetrize=True)
    )
    if pvalue is not None:
        data = data[data['pvalue'] < pvalue]
    data['yagid1'] = id2yagid(data['seqid1'])
    data['yagid2'] = id2yagid(data['seqid2'])

    assert (
        data['yagid1'].str.startswith('YAGID').all() and
        data['yagid2'].str.startswith('YAGID').all()
    )

    ids = ['yagid1', 'yagid2']
    data = data[ids]
    swap = dict(zip(ids, ids[::-1]))
    data = pd.concat([
        data,
        data.rename(columns=swap)
    ])
    data = data[data['yagid1'] < data['yagid2']]

    data = data.drop_duplicates()

    return data
