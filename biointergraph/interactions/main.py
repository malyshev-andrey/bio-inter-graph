import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from tqdm.auto import tqdm


def _fisher_pvalue(n12, n1, n2, n):
    result = fisher_exact(
        [[n12, n1],
         [n2, n]],
        alternative='greater'
    )
    return result.pvalue


def summarize_pairwise(raw_data, symmetrize: bool = False) -> pd.DataFrame:
    if raw_data.shape[1] != 2:
        raise ValueError(f'Two columns are expected, found {raw_data.shape[1]}!')

    id1, id2 = raw_data.columns
    assert (raw_data[id1] != raw_data[id2]).all()
    n_rows = raw_data.shape[0]

    if symmetrize:
        swap = {id1: id2, id2: id1}
        raw_data = pd.concat([
            raw_data,
            raw_data.rename(columns=swap)
        ])

    result = raw_data.groupby([id1, id2], as_index=False).size()
    result['_freq1'] = result.groupby(id1)['size'].transform('sum')
    result['_freq2'] = result.groupby(id2)['size'].transform('sum')

    if symmetrize:
        result = result[result[id1] < result[id2]]

    assert result['size'].sum() == n_rows
    result['_overall'] = result['size'].sum()

    result['PMI'] = np.log2(
        result['size'] * result['_overall']
        / (result['_freq1'] * result['_freq2'])
    )

    result['_freq1'] -= result['size']
    result['_freq2'] -= result['size']
    result['_overall'] -= result['_freq1'] + result['_freq2'] + result['size']

    tqdm.pandas(desc="Fisher's Exact Test calculation")
    result['pvalue'] = result.progress_apply(
        lambda row: _fisher_pvalue(
            row['size'], row['_freq1'],
            row['_freq2'], row['_overall']
        ),
        axis=1
    )
    return result
