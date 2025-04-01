import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, false_discovery_control
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


def summarize_pairwise(
        data: pd.DataFrame,
        ids: list[str], *,
        symmetrize: bool = False,
        fisher_pvalue: bool = True,
        fdr_control: bool = True,
        pmi: bool = True,
        **kwargs
    ) -> pd.DataFrame:
    if len(ids) != 2:
        raise ValueError(f'Two columns are expected, found {ids}!')
    id1, id2 = ids
    assert id1 != id2

    n_pairs = data.shape[0]

    if symmetrize:
        data = pd.concat([
            data,
            data[data[id1] != data[id2]].rename(columns={id1: id2, id2: id1})
        ])

    result = data.groupby([id1, id2], as_index=False).agg(size=(id1, 'size'), **kwargs)
    if fisher_pvalue or pmi:
        result['_freq1'] = result.groupby(id1)['size'].transform('sum')
        result['_freq2'] = result.groupby(id2)['size'].transform('sum')
        result['_overall'] = n_pairs

    if symmetrize:
        result = result[result[id1] <= result[id2]]
    assert result['size'].sum() == n_pairs

    if pmi:
        result['PMI'] = np.log2(
            result['size'] * result['_overall']
            / (result['_freq1'] * result['_freq2'])
        )

    if fisher_pvalue:
        result['_freq1'] -= result['size']
        result['_freq2'] -= result['size']
        result['_overall'] -= result['_freq1'] + result['_freq2'] + result['size']

        tqdm.pandas(desc="Fisher's Exact Test calculation")
        result['pvalue'] = result.progress_apply(_fisher_pvalue, axis=1)
        if fdr_control:
            result['pvalue'] = false_discovery_control(result['pvalue'])

    if fisher_pvalue or pmi:
        result = result.drop(columns=['_freq1', '_freq2', '_overall'])

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
    peaks = peaks.drop_duplicates()

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
