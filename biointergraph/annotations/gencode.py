import re
from ftplib import FTP
from warnings import warn
from typing import Callable

import pandas as pd

from .main import read_feature_table
from ..shared import CHUNKSIZE
from .gff2bed import gff2bed


DOMAIN = 'ftp.ebi.ac.uk'
PATH = 'pub/databases/gencode/Gencode_human'
BLACKLIST = [
    '_README.TXT', '_README_GRCh37_mapping.txt',
    '_README_stats.txt', 'latest_release',
    'stats'
]
DEFAULT_RELEASE = '47'
FORMATS = ('gff3', 'gtf')


def _latest_gencode_release():
    """
    Determines the latest GENCODE release version available on the official FTP server.

    Returns:
        str: The latest GENCODE release version in the format "NN[a-z]?" (e.g., "47", "3b").

    Notes:
        - If the FTP server connection fails, a default version DEFAULT_RELEASE is returned.
        - The function assumes release directories follow the naming pattern `release_NN[a-z]?`.
        - An assertion error is raised if an unexpected naming pattern is encountered.
    """
    try:
        with FTP(DOMAIN) as ftp:
            ftp.login()
            ftp.cwd(PATH)
            releases = [f for f in ftp.nlst() if f not in BLACKLIST]
    except ConnectionRefusedError:
        return DEFAULT_RELEASE

    releases = pd.Series(releases, index=releases)

    regex = r'^release_(?P<version>\d+)(?P<revision>[a-z]?)$'
    assert releases.str.match(regex).all()
    releases = releases.str.extract(regex)

    releases['version'] = releases['version'].astype('int')
    result = releases.sort_values(['version', 'revision']).iloc[-1]

    return result['version'] + result['revision']


def load_gencode_annotation(
        assembly: str, *,
        release: str|int|None = None,
        regions: str, content: str, format: str,
        verbose: bool = True,
        chunksize: int|None = CHUNKSIZE,
        filter_func: Callable[[pd.DataFrame], pd.DataFrame] = lambda df: df,
        **kwargs
    ) -> pd.DataFrame:
    """
    Load GENCODE genome annotations for a specified genome assembly and release.

    This function fetches and parses GENCODE genome annotation files in various formats.
    It supports multiple assemblies and configurations, allowing users to customize
    the regions, content type, file format, and filtering criteria.

    Parameters:
    ----------
    assembly : str
        The genome assembly to use. Valid values are:
        - 'hg19' or 'GRCh37'
        - 'hg38' or 'GRCh38'

    release : str, int, or None, optional (default=None)
        The GENCODE release number (e.g., 47 or '3b'). If not specified, the latest release
        is used.

    regions : str
        Specifies the genomic regions to include. Valid values are:
        - 'chr': Reference chromosomes
        - 'all': Includes chromosomes, patches, haplotypes, and scaffolds
        - 'pri': Primary assembly only

    content : str
        The annotation content type. Valid values are:
        - 'comprehensive': Comprehensive annotation set
        - 'basic': Basic annotation set

    format : str
        File format to fetch. Supported values are gff3, gtf.

    verbose : bool, optional (default=True)
        If True, print the GENCODE URL and feature table shape.

    chunksize : int or None, optional (default=CHUNKSIZE)
        Number of rows per chunk when reading the annotation file.
        If None, the file is read into memory as a single DataFrame.

    filter_func : Callable[[pd.DataFrame], pd.DataFrame], optional (default=lambda df: df)
        A filtering function to apply to the DataFrame after loading.
        Must accept and return a pandas DataFrame.

    **kwargs : dict
        Additional arguments to pass to the `read_feature_table` function.

    Returns:
    -------
    pd.DataFrame
        A DataFrame containing the parsed genome annotation.

    Raises:
    ------
    ValueError
        If invalid values are provided for `assembly`, `release`, `regions`,
        `content`, or `format`.

    Notes:
    -----
    - For `hg19`/`GRCh37`, only 'chr' regions are available. If another value for
        `regions` is specified, it defaults to 'chr' with a warning.
    - The URL for the GENCODE annotation file is constructed dynamically
        based on the input parameters.

    Example:
    -------
    >>> df = load_gencode_annotation(
    ...     assembly='hg38',
    ...     release=39,
    ...     regions='pri',
    ...     content='basic',
    ...     format='gtf',
    ...     verbose=True
    ... )
    GENCODE annotation URL:
        https://.../gencode.v39.primary_assembly.basic.annotation.gtf.gz
    Feature table shape:
        (123456, 10)
    """

    ASSEMBLIES = ('hg19', 'hg38', 'GRCh37', 'GRCh38')
    if assembly not in ASSEMBLIES:
        raise ValueError(
            f'"{assembly}" is not a valid argument. '
            f'Valid arguments are: {", ".join(ASSEMBLIES)}'
        )

    if release is None:
        release = _latest_gencode_release()
    release = str(release)
    if re.fullmatch(r'\d+[a-z]?', release) is None:
        raise ValueError(f'"{release}" is not a valid release name')

    REGIONS = ['chr', 'all', 'pri']
    regions = regions.lower()
    if regions not in REGIONS:
        raise ValueError(
            f'"{regions}" is not a valid argument. '
            f'Valid arguments are: {", ".join(REGIONS)}'
        )

    CONTENTS = ('comprehensive', 'basic')
    content = content.lower()
    if content not in CONTENTS:
        raise ValueError(
            f'"{content}" is not a valid argument. '
            f'Valid arguments are: {", ".join(CONTENTS)}'
        )

    format = format.lower()
    if format not in FORMATS:
        raise ValueError(
            f'"{format}" is not a valid argument. '
            f'Valid arguments are: {", ".join(FORMATS)}'
        )

    is_hg19 = assembly in ('hg19', 'GRCh37')

    if is_hg19 or (regions == 'chr'):
        if regions != 'chr':
            warn(
                'Only CHR (reference chromosomes) regions '
                'are available for hg19 in GENCODE:\n'
                'using regions="CHR" value!'
            )
        regions_suffix = ''
    elif regions == 'all':
        regions_suffix = '.chr_patch_hapl_scaff'
    else:
        assert regions == 'pri'
        regions_suffix = '.primary_assembly'

    full_path = ''.join([
        f'ftp://{DOMAIN}/{PATH}',
        f'/release_{release}',
        '/GRCh37_mapping' if is_hg19 else '',
        f'/gencode.v{release}' + ('lift37' if is_hg19 else ''),
        regions_suffix,
        '.basic' if content == 'basic' else '',
        f'.annotation.{format}.gz'
    ])

    if verbose: print(f'GENCODE annotation URL:\n\t{full_path}')
    result = read_feature_table(full_path, filter_func=filter_func, chunksize=chunksize, **kwargs)
    if verbose: print(f'Feature table shape:\n\t{result.shape}')

    return result


def load_gencode_bed(assembly: str, feature: str) -> pd.DataFrame:
    FEATURES = ('gene', 'transcript')
    if feature not in FEATURES:
        raise ValueError(
            f'"{feature}" is not a valid argument. '
            f'Valid arguments are: {", ".join(FEATURES)}'
        )

    result = []
    for format in FORMATS:
        bed = gff2bed(
            load_gencode_annotation(
                assembly,
                format=format,
                regions='chr' if assembly in ('hg19', 'GRCh37') else 'all',
                content='comprehensive',
                filter_func=lambda df: df[df['type'].eq(feature)]
            ),
            names=feature,
            format=format,
            source='gencode'
        )
        print(f'{format}: {bed.shape[0]}')
        result.append(bed)

    result = pd.concat(result)
    print(f'concat: {result.shape[0]}')
    result = result.drop_duplicates()

    print(f'result: {result.shape[0]}')

    return result
