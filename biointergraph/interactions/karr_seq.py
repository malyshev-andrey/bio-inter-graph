import pandas as pd

from ..shared import memory
from .main import summarize_pairwise
from ..ids_mapping import id2yagid
from .karr_seq_shared import _retrieve_karr_seq_metadata, _load_single_karr_seq


def _load_karr_seq_data(cell_line: str|None = None, **kwargs) -> pd.DataFrame:
    metadata = _retrieve_karr_seq_metadata(cell_line=cell_line)

    result = []
    for _, row in metadata.iterrows():
        data = _load_single_karr_seq(
            row["url"],
            filter_func=lambda df: df[df['seqid1'] != df['seqid2']],
            **kwargs
        )
        for feature in ['dendrimers', 'repl', 'frac', 'cell_line']:
            data[feature] = row[feature]
        result.append(data)
    result = pd.concat(result)

    return result


@memory.cache
def load_karr_seq_data(cell_line: str|None = None, pvalue: float|None):
    result = _load_karr_seq_data(cell_line=cell_line)

    assert (result['seqid1'] != result['seqid2']).all()
    result = result.groupby(
        ['frac', 'cell_line'],
        as_index=False
    ).apply(
        lambda cell_line_frac: summarize_pairwise(
            cell_line_frac,
            ids=['seqid1', 'seqid2'],
            symmetrize=True,
            n_repl=('repl', 'nunique')
        )
    )

    return result

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
