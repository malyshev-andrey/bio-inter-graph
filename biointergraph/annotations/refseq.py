import re
from ftplib import FTP

import pandas as pd

from .main import read_feature_table


def _find_out_latest_refseq_release(
        domain: str = 'ftp.ncbi.nlm.nih.gov',
        path: str = 'genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions',
        assembly: str = 'GRCh38'
    ) -> str:

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
