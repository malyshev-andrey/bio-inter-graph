import re
from ftplib import FTP
from warnings import warn

import pandas as pd

from .main import read_feature_table


def _find_out_latest_gencode_release(
        domain: str = 'ftp.ebi.ac.uk',
        path: str = 'pub/databases/gencode/Gencode_human'
    ) -> str:
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
        verbose: bool = False,
        **kwargs
    ) -> pd.DataFrame:

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

    table = read_feature_table(full_path, **kwargs)
    if verbose:
        print(f'GENCODE annotation URL:\n\t{full_path}')
        print(f'Feature table shape:\n\t{table.shape}')

    return table
