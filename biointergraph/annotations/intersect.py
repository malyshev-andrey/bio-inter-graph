import numpy as np
import pandas as pd
import pyranges as pr

from .ucsc import unify_chr
from ..shared import BED_COLUMNS


BED2RANGES = {
    'chr': 'Chromosome', 'start': 'Start', 'end': 'End',
    'name': 'Name', 'score': 'Score', 'strand': 'Strand'
}
RANGES2BED = {value: key for key, value in BED2RANGES.items()}


def _bed2ranges(bed: pd.DataFrame) -> pr.PyRanges:
    """
    Converts a BED-like pandas DataFrame to a PyRanges object.

    The input DataFrame must have the following columns:
    - `chr` (str): Chromosome name.
    - `start` (int): Start position (0-based).
    - `end` (int): End position (1-based, exclusive).
    - `name` (str): Feature name or identifier.
    - `score` (str): Feature score.
    - `strand` (str): Strand information ('+' or '-').

    Parameters:
        bed (pd.DataFrame): A pandas DataFrame representing genomic intervals in BED format.

    Returns:
        pr.PyRanges: A PyRanges object representing the input intervals.
    """
    result = bed.rename(columns=BED2RANGES)
    result = pr.PyRanges(result)
    return result


def bed_intersect(
        bed1: pd.DataFrame,
        bed2: pd.DataFrame, *,
        strandedness: str|None = 'same',
        unify_chr_assembly: str|None = None,
        jaccard: bool = False,
        overlap: bool = True,
        **kwargs
    ) -> pd.DataFrame:
    """
    Finds intersecting genomic intervals between two BED-like pandas DataFrames.

    The input DataFrames must have the following columns:
    - `chr` (str): Chromosome name.
    - `start` (int): Start position (0-based).
    - `end` (int): End position (1-based, exclusive).
    - `name` (str): Feature name or identifier.
    - `score` (str): Feature score.
    - `strand` (str): Strand information ('+' or '-').

    Intersection is performed using the PyRanges `join` function with strandedness matching ('same' strand only).

    Parameters:
        bed1 (pd.DataFrame): First BED-like DataFrame.
        bed2 (pd.DataFrame): Second BED-like DataFrame.

    Returns:
        pd.DataFrame: A pandas DataFrame containing intersecting intervals, with the following columns:
            - `chr` (str): Chromosome name of the intersecting intervals.
            - `start1`, `end1`, `name1`, `score1`, `strand1`: Columns from the first input DataFrame.
            - `start2`, `end2`, `name2`, `score2`, `strand2`: Columns from the second input DataFrame.
    """

    if unify_chr_assembly is not None:
        bed1['chr'] = unify_chr(bed1['chr'], assembly=unify_chr_assembly)
        bed2['chr'] = unify_chr(bed2['chr'], assembly=unify_chr_assembly)

    default_kwargs = {
        'report_overlap': True,
        'suffix': '_b'
    }
    if strandedness is not None:
        default_kwargs['strandedness'] = strandedness
    default_kwargs.update(kwargs)

    result = _bed2ranges(bed1).join(_bed2ranges(bed2), **default_kwargs).df

    result = result.rename(columns={
        'Chromosome': 'chr',
        'Start': 'start1', 'End': 'end1', 'Name': 'name1',
        'Score': 'score1', 'Strand': 'strand1',
        'Start_b': 'start2', 'End_b': 'end2', 'Name_b': 'name2',
        'Score_b': 'score2', 'Strand_b': 'strand2'
    })

    if jaccard:
        union = (
            np.maximum(result['end1'], result['end2'])
            - np.minimum(result['start1'], result['start2'])
        )
        intersect = (
            np.minimum(result['end1'], result['end2'])
            - np.maximum(result['start1'], result['start2'])
        )
        assert (result['Overlap'] == intersect).all()
        result['jaccard'] = intersect / union

    if not overlap:
        result = result.drop(columns='Overlap')

    return result


def bed_merge(bed: pd.DataFrame, **kwargs) -> pd.DataFrame:
    result = _bed2ranges(bed)

    result = result.merge(**kwargs)

    result = result.df.rename(columns=RANGES2BED)
    return result


def bed_cluster(bed: pd.DataFrame, **kwargs) -> pd.DataFrame:
    result = _bed2ranges(bed)

    result = result.cluster(**kwargs)

    result = result.df.rename(columns=RANGES2BED)
    return result


def best_left_intersect(
        bed1: pd.DataFrame, bed2: pd.DataFrame, *,
        stranded: bool = True,
        unify_chr_assembly: str|None = None,
        jaccard: float|None = None,
        drop_duplicates: bool = True,
        **kwargs
    ) -> pd.DataFrame:

    result = bed_intersect(
        bed1, bed2,
        unify_chr_assembly=unify_chr_assembly,
        strandedness = 'same' if stranded else None,
        jaccard=True, overlap=False,
        apply_strand_suffix=False,
        how='left',
        **kwargs
    )
    result = result.rename(columns={f'{c}1': c for c in BED_COLUMNS})
    assert bed1.columns.isin(result.columns).all()
    left_columns = list(bed1.columns)
    right_columns = [c for c in result.columns if c not in left_columns]

    result.loc[result['start2'].eq(-1), right_columns] = float('nan')

    if drop_duplicates:
        result = result.sort_values('jaccard')
        result = result.drop_duplicates(left_columns, keep='last')
        assert result.shape[0] == bed1.shape[0], f'{result.shape} {bed1.shape} {result["chr"].nunique()} {bed1["chr"].nunique()}'
    else:
        max_jaccard = result.groupby(left_columns, observed=True)['jaccard'].transform('max')
        result = result[(result['jaccard'] == max_jaccard) | result['jaccard'].isna()]

    if jaccard is not None:
        result.loc[result['jaccard'] < jaccard, right_columns] = float('nan')

    return result
