from typing import Callable

import pandas as pd


BED_COLUMNS = [
    'chr', 'start', 'end',
    'name', 'score', 'strand'
]


def _gff2gene_id(ft: pd.DataFrame, *, format: str, source: str) -> pd.Series:
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


def gff2bed(ft: pd.DataFrame, *, names: pd.Series | Callable) -> pd.DataFrame:
    ft = ft.copy()

    try:
        ft['name'] = names(ft)
    except TypeError:
        ft['name'] = names

    return ft[BED_COLUMNS]
