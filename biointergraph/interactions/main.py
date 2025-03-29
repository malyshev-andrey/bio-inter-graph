import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from tqdm.auto import tqdm

from ..annotations import best_left_intersect
from ..ids_mapping import id2yapid, id2yagid


def _fisher_pvalue(row: pd.Series) -> float:
    result = fisher_exact(
        [[row['size'], row['_freq1']],
         [row['_freq2'], row['_overall']]],
        alternative='greater'
    )
    result = result.pvalue
    return result


def summarize_pairwise(pairs: pd.DataFrame, *, symmetrize: bool = False) -> pd.DataFrame:
    if pairs.shape[1] != 2:
        raise ValueError(f'Two columns are expected, found {raw_data.shape[1]}!')
    id1, id2 = pairs.columns
    assert id1 != id2

    n_pairs = pairs.shape[0]

    if symmetrize:
        pairs = pd.concat([
            pairs,
            pairs[pairs[id1] != pairs[id2]].rename(columns={id1: id2, id2: id1})
        ])

    result = pairs.groupby([id1, id2], as_index=False).size()
    result['_freq1'] = result.groupby(id1)['size'].transform('sum')
    result['_freq2'] = result.groupby(id2)['size'].transform('sum')

    if symmetrize:
        result = result[result[id1] <= result[id2]]

    assert result['size'].sum() == n_pairs
    result['_overall'] = result['size'].sum()

    result['PMI'] = np.log2(
        result['size'] * result['_overall']
        / (result['_freq1'] * result['_freq2'])
    )

    result['_freq1'] -= result['size']
    result['_freq2'] -= result['size']
    result['_overall'] -= result['_freq1'] + result['_freq2'] + result['size']

    tqdm.pandas(desc="Fisher's Exact Test calculation")
    result['pvalue'] = result.progress_apply(_fisher_pvalue, axis=1)
    return result


def _annotate_peaks(
        peaks: pd.DataFrame,
        annotation: pd.DataFrame, *,
        assembly: str,
        desc: str|None = None,
        stranded: bool = True,
        convert_ids: bool = False,
        **kwargs
    ) -> pd.DataFrame:
    result = best_left_intersect(
        peaks, annotation,
        stranded=stranded,
        unify_chr_assembly=assembly,
        **kwargs
    )

    no_intersect = result['jaccard'].isna()
    if desc is not None:
        print(f'{desc} peaks without intersections: {no_intersect.sum()}')
    result = result[~no_intersect]

    result = result[['name', 'name2']]
    result.columns = 'source', 'target'

    if convert_ids:
        result['source'] = id2yapid('SYMBOL:' + result['source'], strict=True)
        result['target'] = id2yagid(result['target'], strict=True)

        result = result.drop_duplicates()

    return result
