from ftplib import FTP
from typing import Callable

import pandas as pd

from .main import read_feature_table
from .gff2bed import gff2bed
from ..shared import CHUNKSIZE, memory


DOMAIN = 'ftp.ncbi.nlm.nih.gov'
PATH = 'genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions'
FORMATS = ('gtf', 'gff')

def _latest_refseq_release(assembly: str) -> str:
    """
    Determines the latest RefSeq release version available for a specified genome assembly.

    This function connects to the NCBI FTP server to list all available RefSeq assembly versions
    for Homo sapiens. It filters and extracts the relevant version number for the specified assembly
    (e.g., "GRCh38" or "T2T") and returns the identifier of the latest version.

    Parameters
    ----------
    assembly : str
        The name of the genome assembly for which to determine the latest RefSeq release version.
        Valid values are extracted dynamically from available assemblies on the FTP server.

    Returns
    -------
    str
        The identifier of the latest RefSeq release version for the specified genome assembly.

    Raises
    ------
    ValueError
        If the specified assembly is not found among the available assemblies on the server.
    AssertionError
        If the server data does not match the expected format or contains inconsistencies.

    Notes
    -----
    - The function expects the FTP directory structure and file naming conventions to match
      those used by NCBI RefSeq.

    Examples
    --------
    >>> _latest_refseq_release("GRCh38")
    'GCF_000001405.40_GRCh38.p14'
    """
    with FTP(DOMAIN) as ftp:
        ftp.login()
        ftp.cwd(PATH)
        releases = ftp.nlst()

    releases.remove('suppressed')
    releases = pd.Series(releases, index=releases)

    regex = r'^GCF_\d{9}\.(?P<version>\d+)_(?P<assembly>NCBI|GRCh37|GRCh38|T2T)'
    assert releases.str.match(regex).all()
    releases = releases.str.extract(regex)

    assemblies = releases['assembly'].unique()
    if assembly not in assemblies:
        raise ValueError(
            f'"{assembly}" is not a valid argument. '
            f'Valid arguments are: {", ".join(assemblies)}'
        )

    releases = releases.loc[releases['assembly'] == assembly, 'version']
    assert releases.is_unique and releases.str.isdigit().all()
    releases = releases.astype('int')

    result = releases.idxmax()

    return result


def load_refseq_annotation(
        assembly: str, *,
        format: str,
        verbose: bool = True,
        chunksize: int | None = CHUNKSIZE,
        filter_func: Callable[[pd.DataFrame], pd.DataFrame] = lambda df: df,
        **kwargs
    ) -> pd.DataFrame:
    """
    Load RefSeq genome annotations for a specified genome assembly and format.

    Parameters
    ----------
    assembly : str
        The genome assembly to load annotations for. Valid options are:
        - 'hg38' or 'GRCh38' (GRCh38 assembly)
        - 'hg19' or 'GRCh37' (GRCh37 assembly)
        - 'T2T' (T2T assembly)
    format : str
        The format of the annotation file. Valid options are:
        - 'gtf'
        - 'gff'
    verbose : bool, optional
        If True, print progress messages and details about the loaded data. Defaults to True.
    chunksize : int or None, optional
        Number of rows to read in each chunk. If None, the entire file is loaded into memory. Defaults to CHUNKSIZE.
    filter_func : Callable[[pd.DataFrame], pd.DataFrame], optional
        A callable that filters or processes the loaded DataFrame. It takes a DataFrame as input
        and returns a processed DataFrame. Defaults to an identity function (no filtering).
    **kwargs : dict, optional
        Additional arguments to pass to the internal data loading function (`read_feature_table`).

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the RefSeq genome annotation data.

    Raises
    ------
    ValueError
        If an invalid `assembly` or `format` is specified.

    Notes
    -----
    - This function uses the `_latest_refseq_release` function to determine the latest RefSeq release for the
        specified assembly.
    - The annotation data is downloaded as a compressed file from a predefined URL structure.
    - The `read_feature_table` function is responsible for parsing the downloaded data.

    Examples
    --------
    Load GRCh38 GTF annotation with verbose output:

    >>> df = load_refseq_annotation(assembly='GRCh38', format='gtf', verbose=True)

    Load hg19 GFF annotation, process chunks of 10000 rows, and filter genes:

    >>> df = load_refseq_annotation(
    ...     assembly='hg19',
    ...     format='gff',
    ...     chunksize=1e4,
    ...     filter_func=lambda df: df[df['type'] == 'gene']
    ... )
    """
    ASSEMBLIES = {
        'hg38': 'GRCh38', 'GRCh38': 'GRCh38',
        'hg19': 'GRCh37', 'GRCh37': 'GRCh37',
        'T2T': 'T2T'
    }
    if assembly not in ASSEMBLIES:
        raise ValueError(
            f'"{assembly}" is not a valid argument. '
            f'Valid arguments are: {", ".join(ASSEMBLIES)}'
        )

    format = format.lower()
    if format not in FORMATS:
        raise ValueError(
            f'"{format}" is not a valid argument. '
            f'Valid arguments are: {", ".join(FORMATS)}'
        )

    release = _latest_refseq_release(ASSEMBLIES[assembly])
    full_path = '/'.join([f'ftp://{DOMAIN}', PATH, f'{release}/{release}_genomic.{format}.gz'])

    if verbose: print(f'RefSeq annotation URL:\n\t{full_path}')
    result = read_feature_table(full_path, chunksize=chunksize, filter_func=filter_func, **kwargs)
    if verbose: print(f'Feature table shape:\n\t{result.shape}')

    return result


@memory.cache
def load_refseq_bed(assembly: str, feature: str) -> pd.DataFrame:
    FEATURES = ('gene', 'transcript')
    if feature not in FEATURES:
        raise ValueError(
            f'"{feature}" is not a valid argument. '
            f'Valid arguments are: {", ".join(FEATURES)}'
        )

    shapes = []
    result = []
    for format in FORMATS:
        bed = gff2bed(
            load_refseq_annotation(
                assembly,
                format=format,
                filter_func=lambda df: df[df['type'].eq(feature)]
            ),
            names=feature,
            format=format,
            source='refseq'
        )
        shapes.append(f'{format}: {bed.shape[0]}')
        result.append(bed)

    shapes = [' | '.join(shapes)]

    result = pd.concat(result)
    shapes.append(str(result.shape[0]))

    result = result.drop_duplicates()
    shapes.append(str(result.shape[0]))

    print(f'RefSeq {feature}s to BED-format: {" -> ".join(shapes)}')
    return result
