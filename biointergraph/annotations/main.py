from pathlib import Path
from warnings import warn
from typing import IO, Callable

import pandas as pd
from tqdm.auto import tqdm

from ..shared import GFF_COLUMNS


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
        filter_func: Callable[[pd.DataFrame], pd.DataFrame] = lambda df: df,
        chunksize: int | None = None,
        validation: bool = True,
        **kwargs
    ) -> pd.DataFrame:
    """
    Reads a GFF, GTF, or GFF3 feature table from a file, URL, or other sources,
    and applies optional filtering and validation.

    This function loads a feature table, optionally in chunks, applies a user-defined filtering
    function, and validates the resulting table. By default, the feature table is expected to
    be tab-delimited with specific columns defined in `GFF_COLUMNS`.

    Args:
        filepath_or_buffer (str | Path | IO[str]): The input source of the feature table.
            This can be:
            - A local file path (str or Path).
            - A URL pointing to a remote file (str).
            - A file-like object (IO[str]).
        filter_func (Callable[[pd.DataFrame], pd.DataFrame], optional): A function to filter
            or transform the feature table after loading. Defaults to the identity function.
        chunksize (int | None, optional): The number of rows to read per chunk. If None,
            the entire file is read at once. Defaults to None.
        validation (bool, optional): Whether to validate the feature table after reading.
            Defaults to True.
        **kwargs: Additional arguments passed to `pandas.read_csv`, such as `sep` or `dtype`.

    Returns:
        pd.DataFrame: The processed feature table with valid 'start' and 'end' columns
        converted to integers.

    Raises:
        ValueError: If the file cannot be read due to formatting issues or missing data.

    Notes:
        - By default, the file is read as a tab-delimited text file with `GFF_COLUMNS` as the header.
        - The `filepath_or_buffer` argument supports any source readable by `pandas.read_csv`,
            including URLs and file-like objects.
        - If `chunksize` is specified, the file is read and processed in chunks, improving
            performance for large files.

    Examples:
        # Reading a GFF file from a local path with validation
        df = read_feature_table('features.gff')

        # Reading a GFF file from a URL
        df = read_feature_table('https://example.com/features.gff')

        # Reading with custom filtering
        df = read_feature_table('features.gff', filter_func=lambda df: df[df['type'] == 'gene'])

        # Reading in chunks
        df = read_feature_table('features.gff', chunksize=1000)
    """

    read_csv_kwargs = dict(
        sep='\t', comment='#',
        header=None, names=GFF_COLUMNS,
        dtype='str'
    )
    read_csv_kwargs.update(kwargs)

    if chunksize is None:
        result = filter_func(pd.read_csv(filepath_or_buffer, **read_csv_kwargs))
    else:
        result = []

        if not isinstance(filepath_or_buffer, str):
            desc = 'Reading feature table rows: '
        elif len(filepath_or_buffer) <= 60:
            desc = filepath_or_buffer
        else:
            desc = filepath_or_buffer[30] + ' ... ' + filepath_or_buffer[-30:]

        with tqdm(desc=desc) as progress_bar:
            for chunk in pd.read_csv(filepath_or_buffer, chunksize=chunksize, **read_csv_kwargs):
                progress_bar.update(chunk.shape[0])
                result.append(filter_func(chunk))
        result = pd.concat(result)

    if validation:
        result = _validate_feature_table(result)

    result['start'] = result['start'].astype('int')
    result['end'] = result['end'].astype('int')

    return result
