from ftplib import FTP
import gzip
from io import BytesIO

import requests
import pandas as pd

from ..annotations import unify_chr

DOMAIN = 'ftp.ensembl.org'


def _latest_ensembl_release() -> str:
    with FTP(DOMAIN) as ftp:
        ftp.login()
        ftp.cwd('pub')
        result = ftp.nlst()

    result = pd.Series(result)
    release_regex = r'^release-(?P<version>\d+)$'
    result = result[result.str.match(release_regex)]
    result = result.str.extract(release_regex)['version']

    result.index = result
    assert result.str.isdigit().all()
    result = result.astype('int')
    assert result.max() >= 113
    result = result.idxmax()

    return result


def _ensembl_mysql_prefix(release: str|None = None) -> str:
    if release is None:
        release = _latest_ensembl_release()

    with FTP(DOMAIN) as ftp:
        ftp.login()
        ftp.cwd(f'pub/release-{release}/mysql')
        versions = ftp.nlst()

    versions = pd.Series(versions)
    regex = r'^homo_sapiens_core_(?P<release>\d+)_(?P<version>\d+)(?P<subversion>[a-z]*)$'
    versions = versions[versions.str.match(regex)]
    versions = versions.str.extract(regex)

    assert versions['release'].eq(release).all()

    versions['version'] = versions['version'].astype('int')
    versions = versions.sort_values(['version', 'subversion']).iloc[-1]
    result = f'homo_sapiens_core_{release}_{versions["version"]}{versions["subversion"]}'

    url = f'https://{DOMAIN}/pub/release-{release}/mysql/{result}'
    assert requests.get(url).status_code == 200

    return result


def _retrieve_ensembl_schema(table, *, release: str|None = None) -> list[str]:
    prefix = _ensembl_mysql_prefix(release)
    path = f'pub/current_mysql/{prefix}/{prefix}.sql.gz'

    with FTP(DOMAIN) as ftp, BytesIO() as stream:
        ftp.login()
        ftp.retrbinary(f"RETR {path}", stream.write)
        stream.seek(0)

        with gzip.GzipFile(fileobj=stream) as input_file:
            name = None
            name2schema = {}
            for line in input_file:
                line = line.decode()

                if line.startswith('CREATE TABLE'):
                    name = line.split('`')[1]
                    assert name not in name2schema
                    name2schema[name] = []

                elif line.startswith('  `'):
                    assert name is not None
                    column = line.split('`')[1]
                    name2schema[name].append(column)

    assert table in name2schema
    result = name2schema[table]

    return result


def fetch_ensembl_table(table: str, *, release: str|None = None) -> pd.DataFrame:
    prefix = _ensembl_mysql_prefix(release)
    url = f'ftp://ftp.ensembl.org/pub/current_mysql/{prefix}/{table}.txt.gz'
    result = pd.read_csv(url, sep='\t', header=None, dtype='str')

    columns = _retrieve_ensembl_schema(table, release=release)
    assert result.shape[1] == len(columns)
    result.columns = columns

    return result


def ensembl_transcript_id_info(release: str|None = None) -> pd.DataFrame:
    result = fetch_ensembl_table('transcript', release=release)
    chr_names = fetch_ensembl_table('seq_region', release=release)
    chr_names = chr_names[['seq_region_id', 'name']]

    result = result.merge(chr_names, how='left', validate='many_to_one')

    columns = {
        'seq_region_start': 'start',
        'seq_region_end': 'end',
        'seq_region_strand': 'strand',
        'biotype': 'biotype',
        'stable_id': 'ensembl_transcript_id',
        'name': 'chr'
    }
    result = result.rename(columns=columns)
    result = result[columns.values()]

    result['chr'] = unify_chr(result['chr'], assembly='hg38')
    assert result['chr'].eq('chrM').any() and not result['chr'].eq('chrMT').any()

    return result


def ensembl_gene_id_info(release: str|None = None) -> pd.DataFrame:
    result = fetch_ensembl_table('gene', release=release)
    chr_names = fetch_ensembl_table('seq_region', release=release)
    chr_names = chr_names[['seq_region_id', 'name']]

    result = result.merge(chr_names, how='left', validate='many_to_one')
    columns = {
        'seq_region_start': 'start',
        'seq_region_end': 'end',
        'seq_region_strand': 'strand',
        'biotype': 'biotype',
        'stable_id': 'ensembl_gene_id',
        'name': 'chr'
    }
    result = result.rename(columns=columns)
    result = result[columns.values()]

    result['chr'] = unify_chr(result['chr'], assembly='hg38')
    assert result['chr'].eq('chrM').any() and not result['chr'].eq('chrMT').any()

    return result
