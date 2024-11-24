import re
from ftplib import FTP

import pandas as pd

from .main import read_feature_table


def _find_out_latest_refseq_release(
        domain: str = 'ftp.ncbi.nlm.nih.gov',
        path: str = 'genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions',
        assembly: str = 'GRCh38'
    ) -> str:
    """
    Determines the latest RefSeq release version available for a specified assembly.

    Args:
        domain (str): FTP server domain to connect to. Defaults to 'ftp.ncbi.nlm.nih.gov'.
        path (str): Path on the FTP server where RefSeq assemblies are stored.
            Defaults to 'genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions'.
        assembly (str): Genome assembly to check for. Valid options are:
            - 'GRCh37': Genome Reference Consortium Human Build 37.
            - 'GRCh38': Genome Reference Consortium Human Build 38.
            - 'T2T': Telomere-to-Telomere assembly CHM13vX.
            Defaults to 'GRCh38'.

    Returns:
        str: The latest RefSeq release directory name corresponding to the specified assembly.

    Raises:
        AssertionError: If the `assembly` argument is not valid.

    Notes:
        - The function identifies release versions using assembly-specific patterns.
        - Release directories are sorted based on version numbers to select the latest one.
        - The return value includes the directory name for the latest release.
    """


    ftp = FTP(domain)
    ftp.login()
    ftp.cwd(path)

    search_pattern = {
        'GRCh37': r'^GCF_000001405\.(\d+)_GRCh37(?:\.p(\d+))?$',
        'GRCh38': r'^GCF_000001405\.(\d+)_GRCh38(?:\.p(\d+))?$',
        'T2T': r'^GCF_009914755\.(\d+)_T2T-CHM13v(\d+)\.(\d+)$'
    }

    assert assembly in search_pattern

    releases = []

    for name in ftp.nlst():
        if assembly in name:
            version = re.search(search_pattern[assembly], name).groups()
            version = filter(lambda s: s is not None, version)
            version = map(int, version)
            releases.append((*version, name))

    return sorted(releases)[-1][-1]


def load_refseq_annotation(
        assembly: str = 'GRCh38', *,
        format: str = 'gff',
        verbose: bool = False,
        domain='ftp.ncbi.nlm.nih.gov',
        path='genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions',
        **kwargs
    ) -> pd.DataFrame:
    """
    Loads RefSeq genome annotations for a specified assembly.

    Args:
        assembly (str): Genome assembly to use. Valid options are:
            - 'hg38' or 'GRCh38': Genome Reference Consortium Human Build 38.
            - 'hg19' or 'GRCh37': Genome Reference Consortium Human Build 37.
            - 'T2T': Telomere-to-Telomere assembly CHM13vX.
            Defaults to 'GRCh38'.
        format (str): File format for the annotation. Valid options are:
            - 'gff': General Feature Format (default).
            - 'gtf': Gene Transfer Format.
        verbose (bool): If True, prints the RefSeq URL and the feature table shape.
            Defaults to False.
        domain (str): FTP server domain to connect to. Defaults to 'ftp.ncbi.nlm.nih.gov'.
        path (str): Path on the FTP server where RefSeq assemblies are stored.
            Defaults to 'genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions'.
        **kwargs: Additional keyword arguments passed to `read_feature_table`.

    Returns:
        pd.DataFrame: A feature table containing RefSeq annotation data.

    Raises:
        ValueError: If an invalid `assembly` argument is provided.

    Notes:
        - The function determines the latest release version for the specified assembly
            using `_find_out_latest_refseq_release`.
        - The annotation data is fetched from the dynamically constructed URL based on
            the selected assembly and format.
        - If `verbose=True`, the constructed URL and feature table dimensions are printed.
    """

    ASSEMBLIES = {
        'hg38': 'GRCh38',
        'hg19': 'GRCh37',
        'GRCh37': 'GRCh37',
        'GRCh38': 'GRCh38',
        'T2T': 'T2T'
    }

    if assembly not in ASSEMBLIES:
        raise ValueError(f'"{assembly}" is not a valid argument. Valid arguments are: {", ".join(ASSEMBLIES)}')

    v = _find_out_latest_refseq_release(assembly=ASSEMBLIES[assembly], domain=domain, path=path)

    full_path = '/'.join([f'https://{domain}', path, f'{v}/{v}_genomic.{format}.gz'])

    table = read_feature_table(full_path, **kwargs)

    if verbose:
        print(f'RefSeq annotation URL:\n\t{full_path}')
        print(f'Feature table shape:\n\t{table.shape}')

    return table
