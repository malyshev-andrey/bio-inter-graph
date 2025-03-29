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
def load_karr_seq_data(cell_line: str|None, pvalue: float|None = None) -> pd.DataFrame:
    result = _load_karr_seq_data(cell_line=cell_line)

    assert (result['seqid1'] != result['seqid2']).all()
    result = result.groupby(
        ['cell_line', 'frac']
    ).apply(
        lambda cell_line_frac: summarize_pairwise(
            cell_line_frac,
            ids=['seqid1', 'seqid2'],
            symmetrize=True,
            pmi=False
        )
    )

    if pvalue is not None:
        result = result[result['pvalue'] < pvalue]

    result = pd.DataFrame({
        'yagid1': id2yagid(result['seqid1'], strict=True),
        'yagid2': id2yagid(result['seqid2'], strict=True)
    })

    swap_mask = result['yagid1'] > result['yagid2']
    result.loc[swap_mask, ['yagid1', 'yagid2']] = result.loc[swap_mask, ['yagid2', 'yagid1']].values
    result = result[result['yagid1'] < result['yagid2']]

    result = result.drop_duplicates()

    return result
