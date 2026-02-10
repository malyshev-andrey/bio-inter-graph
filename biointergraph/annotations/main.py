from pathlib import Path
from warnings import warn
from typing import IO, Callable

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from ..shared import GFF_COLUMNS, _read_tsv


def _validate_feature_table(ft: pd.DataFrame) -> pd.DataFrame:
    """
    Validates the feature table to ensure required columns have correct formats.

    This function checks that the 'start' and 'end' columns contain numeric values,
    and the 'strand' column contains only valid strand identifiers ('+', '-', or '.').
    Invalid rows are filtered out, and warnings are issued for the first few
    incorrect entries in each column.

    Args:
        ft (pd.DataFrame): A DataFrame representing the feature table, expected to
            have at least 'start', 'end', and 'strand' columns.

    Returns:
        pd.DataFrame: A filtered DataFrame containing only valid rows.

    Warnings:
        - Warns about incorrect 'start' values, showing up to the first 5 invalid entries.
        - Warns about incorrect 'end' values, showing up to the first 5 invalid entries.
        - Warns about incorrect 'strand' values, showing up to the first 5 invalid entries.
    """

    is_valid_start = ft['start'].str.isdigit()
    if not is_valid_start.all():
        incorrect = " ".join(ft['start'][~is_valid_start].head(5))
        warn(
            f'Some start coordinates seem to be incorrect: {incorrect}.\n'
            'Skipping these lines!'
        )

    is_valid_end = ft['end'].str.isdigit()
    if not is_valid_end.all():
        incorrect = " ".join(ft['end'][~is_valid_end].head(5))
        warn(
            f'Some end coordinates seem to be incorrect: {incorrect}.\n'
            'Skipping these lines!'
        )

    is_valid_strand = ft['strand'].isin(['+', '-', '.'])
    if not is_valid_strand.all():
        incorrect = " ".join(ft['strand'][~is_valid_strand].head(5))
        warn(
            f'Some strand values seem to be incorrect: {incorrect}.\n'
            'Skipping these lines!'
        )

    ft = ft[is_valid_strand & is_valid_end & is_valid_start]

    return ft


def read_feature_table(
        filepath_or_buffer: str | Path | IO[str], *,
        validation: bool = True,
        **kwargs
    ) -> pd.DataFrame:

    result = _read_tsv(
        filepath_or_buffer,
        comment='#',
        header=None,
        names=GFF_COLUMNS,
        use_cache=True,
        **kwargs
    )

    if validation:
        result = _validate_feature_table(result)

    assert result['start'].str.isdigit().all()
    result['start'] = result['start'].astype('int')

    assert result['end'].str.isdigit().all()
    result['end'] = result['end'].astype('int')

    assert (result['start'] <= result['end']).all()

    return result


def sanitize_bed(
        bed: pd.DataFrame, *,
        stranded: bool = True,
        inplace: bool = False
    ) -> pd.DataFrame|None:
    if not inplace:
        bed = bed.copy()

    assert bed['start'].str.isdigit().all()
    bed['start'] = bed['start'].astype('int')

    assert bed['end'].str.isdigit().all()
    bed['end'] = bed['end'].astype('int')

    assert (bed['start'] < bed['end']).all()

    if stranded and 'strand' in bed.columns:
        strand_values = {'+', '-'} if stranded else {'+', '-', '.'}
        assert bed['strand'].isin(strand_values).all()

    return None if inplace else bed


def _row2intervals_factory(bin_size: int) -> Callable[[pd.Series], tuple[np.array, np.array]]:
    def _row2intervals(row: pd.Series) -> tuple[np.array, np.array]:
        starts = np.linspace(
            row['start'], row['end'],
            num=1 + max(1, (row['end'] - row['start']) // bin_size),
            dtype='int'
        )

        return starts[:-1], np.roll(starts, -1)[:-1]

    return _row2intervals


def _split_annotation_into_bins(annotation: pd.DataFrame, bin_size: int) -> pd.DataFrame:
    missing_columns = [c for c in ('start', 'end') if c not in annotation.columns]
    if missing_columns:
        raise ValueError(f'Missing required columns: {", ".join(missing_columns)}')
    assert (annotation['end'] > annotation['start']).all()

    result = annotation.copy()

    tqdm.pandas(desc='Splitting annotation intervals')
    result[['start', 'end']] = result.progress_apply(
        _row2intervals_factory(bin_size),
        result_type='expand',
        axis=1
    )
    result = result.explode(['start', 'end'])
    result = result.astype({'start': 'int', 'end': 'int'})

    length = result['end'] - result['start']
    assert length.between(1, 2 * bin_size - 1).all()
    assert length.sum() == (annotation['end'] - annotation['start']).sum()

    return result
