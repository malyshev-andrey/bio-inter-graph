import re
from ftplib import FTP
from warnings import warn

import pandas as pd

from .main import read_feature_table


def _find_out_latest_gencode_release(
        domain: str = 'ftp.ebi.ac.uk',
        path: str = 'pub/databases/gencode/Gencode_human'
    ) -> str:
    """
    Determines the latest GENCODE release version available on the official FTP server.

    Args:
        domain (str): FTP server domain to connect to. Defaults to 'ftp.ebi.ac.uk'.
        path (str): Path on the FTP server where GENCODE releases are stored.
            Defaults to 'pub/databases/gencode/Gencode_human'.

    Returns:
        str: The latest GENCODE release version in the format "NN[a-z]" (e.g., "47", "47a").

    Notes:
        - If the FTP server connection fails, a default version "47" is returned.
        - The function assumes release directories follow the naming pattern `release_NN[a-z]?`.
        - An assertion error is raised if an unexpected naming pattern is encountered.
    """

    try:
        ftp = FTP(domain)
    except ConnectionRefusedError:
        return '47'

    ftp.login()
    ftp.cwd(path)
    releases = []
    for name in ftp.nlst():
        if name.startswith('release_'):
            is_match = re.fullmatch(r'^release_\d+[a-z]?$', name) is not None
            assert is_match, f'Unexpected GENCODE release name: {name}'
            result = re.search(r'^release_(\d+)([a-z]?)$', name)
            version = int(result.group(1))
            revision = result.group(2)
            releases.append((version, revision))
    version, revision = sorted(releases)[-1]
    result = str(version) + revision
    return result


def load_gencode_annotation(
        *, assembly: str = 'hg38',
        release: str|int|None = None,
        regions: str = 'CHR',
        content: str = 'comprehensive',
        format: str = 'gff3',
        verbose: bool = True,
        **kwargs
    ) -> pd.DataFrame:
    """
    Loads GENCODE genome annotations for a specified assembly and release.

    Args:
        assembly (str): Genome assembly to use. Valid options are 'hg19', 'hg38',
            'GRCh37', and 'GRCh38'. Defaults to 'hg38'.
        release (str | int | None): GENCODE release version (e.g., '47' or '47a').
            If None, the latest release is determined automatically. Defaults to None.
        regions (str): Specifies the region of the genome to download. Valid options are:
            - 'chr': Reference chromosomes only.
            - 'all': Includes additional patches, haplotypes, and scaffolds.
            - 'pri': Primary assembly only.
            Defaults to 'CHR'.
        content (str): Annotation content type. Valid options are:
            - 'comprehensive': Includes all annotated features.
            - 'basic': Limited to core transcripts.
            Defaults to 'comprehensive'.
        format (str): File format for the annotation. Valid options are 'gff3' and 'gtf'.
            Defaults to 'gff3'.
        verbose (bool): If True, prints the GENCODE URL and the feature table shape.
            Defaults to True.
        **kwargs: Additional keyword arguments passed to `read_feature_table`.

    Returns:
        pd.DataFrame: A feature table containing GENCODE annotation data.

    Raises:
        ValueError: If invalid arguments are provided for `assembly`, `release`,
            `regions`, `content`, or `format`.

    Notes:
        - For 'hg19' or 'GRCh37', only CHR (reference chromosomes) regions are supported.
        - The function constructs the GENCODE annotation URL dynamically based on
            the inputs and downloads the data using the `read_feature_table` function.
        - The full URL and feature table shape are printed if `verbose=True`.
    """

    ASSEMBLIES = ('hg19', 'hg38', 'GRCh37', 'GRCh38')
    if assembly not in ASSEMBLIES:
        raise ValueError(f'"{assembly}" is not a valid argument. Valid arguments are: {", ".join(ASSEMBLIES)}')

    if release is None:
        release = _find_out_latest_gencode_release()
    release = str(release)
    if re.fullmatch(r'\d+[a-z]?', release) is None:
        raise ValueError(f'"{release}" is not a valid release name')

    REGIONS = ['chr', 'all', 'pri']
    regions = regions.lower()
    if regions not in REGIONS:
        raise ValueError(f'"{regions}" is not a valid argument. Valid arguments are: {", ".join(REGIONS)}')

    CONTENTS = ('comprehensive', 'basic')
    content = content.lower()
    if content not in CONTENTS:
        raise ValueError(f'"{content}" is not a valid argument. Valid arguments are: {", ".join(CONTENTS)}')

    FORMATS = ('gff3', 'gtf')
    format = format.lower()
    if format not in FORMATS:
        raise ValueError(f'"{format}" is not a valid argument. Valid arguments are: {", ".join(FORMATS)}')

    human_path = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human'

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
        human_path,
        f'/release_{release}',
        '/GRCh37_mapping' if is_hg19 else '',
        f'/gencode.v{release}' + ('lift37' if is_hg19 else ''),
        regions_suffix,
        '.basic' if content == 'basic' else '',
        f'.annotation.{format}.gz'
    ])

    if verbose: print(f'GENCODE annotation URL:\n\t{full_path}')
    table = read_feature_table(full_path, **kwargs)
    if verbose: print(f'Feature table shape:\n\t{table.shape}')

    return table
