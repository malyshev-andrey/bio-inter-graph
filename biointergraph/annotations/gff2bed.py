from typing import Callable

import numpy as np
import pandas as pd

from ..shared import BED_COLUMNS


def _gff2gene_id(ft: pd.DataFrame, *, format: str, source: str) -> pd.Series:
    """
    Extract gene IDs from the 'attributes' column of a feature table in GFF-like format.

    Parameters:
        ft (pd.DataFrame): Input feature table with columns matching `GFF_COLUMNS`.
        format (str): Format of the input data ('gff3', 'gtf', or 'gff').
        source (str): Data source ('refseq' or 'gencode').

    Returns:
        pd.Series: Extracted gene IDs as a Pandas Series.

    Raises:
        ValueError: If `format` or `source` is invalid or unsupported for the given source.
        AssertionError: If the 'attributes' column contains multiple matches for gene IDs
                        or if 'gene' rows have missing gene IDs.

    Notes:
        - The function uses regex patterns tailored to the specified `source` and `format`.
        - Only GFF, GTF, and GFF3 formats are supported.
        - Input feature table is expected to have 9 columns: ['chr', 'source', 'type', 'start',
            'end', 'score', 'strand', 'phase', 'attributes'].
    """

    FORMATS = ('gff3', 'gtf', 'gff')
    format = format.lower()
    if format not in FORMATS:
        raise ValueError(
            f'"{format}" is not a valid argument. '
            f'Valid arguments are: {", ".join(FORMATS)}'
        )

    SOURCES = ('refseq', 'gencode')
    source = source.lower()
    if source not in SOURCES:
        raise ValueError(
            f'"{source}" is not a valid argument. '
            f'Valid arguments are: {", ".join(SOURCES)}'
        )

    regex = None
    if source == 'refseq':
        if format == 'gtf':
            regex = r'db_xref "GeneID:(\d+)"'
        elif format == 'gff':
            regex = r'GeneID:(\d+)'
        else:
            raise ValueError('Only gff and gtf formats are expected for refseq source!')
    elif source == 'gencode':
        if format == 'gff3':
            regex = r'gene_id=([^;]+)(?:;|$)'
        elif format == 'gtf':
            regex = r'gene_id "([^"]+)"'
        else:
            raise ValueError('Only gff3 and gtf formats are expected for gencode source!')
    else:
        assert source in SOURCES

    assert (ft['attributes'].str.count(regex) <= 1).all()

    result = ft['attributes'].str.extract(regex)[0]
    assert not result[ft['type'].eq('gene')].isna().any()

    return result


def _gff2transcript_id(ft: pd.DataFrame, *, format: str, source: str) -> pd.Series:
    """
    Extract transcript IDs from the 'attributes' column of a feature table in GFF-like format.

    Parameters:
        ft (pd.DataFrame): Input feature table with columns matching `GFF_COLUMNS`.
        format (str): Format of the input data ('gff3', 'gtf', or 'gff').
        source (str): Data source ('refseq' or 'gencode').

    Returns:
        pd.Series: Extracted transcript IDs as a Pandas Series.

    Raises:
        ValueError: If `format` or `source` is invalid or unsupported for the given source.
        AssertionError: If the 'attributes' column contains multiple matches for transcript IDs
                        or if 'transcript' rows have missing transcript IDs.

    Notes:
        - The function uses regex patterns tailored to the specified `source` and `format`.
        - Only GFF, GTF, and GFF3 formats are supported.
        - Input feature table is expected to have 9 columns: ['chr', 'source', 'type', 'start',
            'end', 'score', 'strand', 'phase', 'attributes'].
    """

    FORMATS = ('gff3', 'gtf', 'gff')
    format = format.lower()
    if format not in FORMATS:
        raise ValueError(
            f'"{format}" is not a valid argument. '
            f'Valid arguments are: {", ".join(FORMATS)}'
        )

    SOURCES = ('refseq', 'gencode')
    source = source.lower()
    if source not in SOURCES:
        raise ValueError(
            f'"{source}" is not a valid argument. '
            f'Valid arguments are: {", ".join(SOURCES)}'
        )

    regex = None
    if source == 'refseq':
        if format == 'gtf':
            regex = r'transcript_id "([^"]+)"'
        elif format == 'gff':
            regex = r'transcript_id=([^;]+)(?:;|$)'
        else:
            raise ValueError('Only gff and gtf formats are expected for refseq source!')
    elif source == 'gencode':
        if format == 'gff3':
            regex = r'transcript_id=([^;]+)(?:;|$)'
        elif format == 'gtf':
            regex = r'transcript_id "([^"]+)"'
        else:
            raise ValueError('Only gff3 and gtf formats are expected for gencode source!')
    else:
        assert source in SOURCES

    assert (ft['attributes'].str.count(regex) <= 1).all()

    result = ft['attributes'].str.extract(regex)[0]
    assert not result[ft['type'].eq('transcript')].isna().any()

    return result


def gff2bed(
        ft: pd.DataFrame, *,
        names: pd.Series | Callable[[pd.DataFrame], pd.Series] | str | None = None,
        **kwargs
    ) -> pd.DataFrame:
    """
    Convert a feature table in GFF-like format to a BED-like DataFrame.

    Parameters:
        ft (pd.DataFrame): Input feature table with columns matching `GFF_COLUMNS`.
        names (pd.Series | Callable[[pd.DataFrame], pd.Series]): Specifies the 'name' column
            in the BED output. If a callable is provided, it should accept the feature table
            and additional keyword arguments, returning a Pandas Series of names.
        **kwargs: Additional arguments passed to the callable function in `names`, if applicable.

    Returns:
        pd.DataFrame: BED-like DataFrame with columns matching `BED_COLUMNS`.

    Notes:
        - The resulting DataFrame has columns: ['chr', 'start', 'end', 'name', 'score', 'strand'].
        - Coordinates are converted to zero-based, half-open intervals (BED format).
        - The 'name' column is derived from the `names` parameter, which can be a Pandas Series
            or a callable function that generates names dynamically.
        - Input feature table is expected to have 9 columns: ['chr', 'source', 'type', 'start',
            'end', 'score', 'strand', 'phase', 'attributes'].
    """

    ft = ft.copy()
    ft['name'] = '.'

    alias2func = {
        'gene': _gff2gene_id,
        'transcript': _gff2transcript_id
    }

    if names is None:
        for alias in alias2func:
            ft['name'] = np.where(
                ft['type'].eq(alias),
                alias2func[alias](ft, **kwargs),
                ft['name']
            )
        invalid_frac = ft['name'].eq('.').mean()
        if invalid_frac > 0:
            print(f'gff2bed: output bed invalid names frac: {invalid_frac:.04f}')
    else:
        if isinstance(names, str):
            if names not in alias2func:
                raise ValueError(
                    f'"{names}" is not a valid argument. '
                    f'Valid arguments are: {", ".join(alias2func)}'
                )
            names = alias2func[names]
        ft['name'] = names(ft, **kwargs) if callable(names) else names

    ft['start'] = ft['start'].astype('int') - 1
    ft['end'] = ft['end'].astype('int')

    assert (ft['start'] < ft['end']).all()

    return ft[BED_COLUMNS]
