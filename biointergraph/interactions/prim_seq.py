from io import StringIO

import requests
import pandas as pd

from ..shared import _read_tsv, memory
from ..ids_mapping import id2yagid, id2yapid


LINK1 = 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8332nnn/GSM8332740/suppl/GSM8332740_K562_1_chimericReads.csv.gz'
LINK2 = 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8332nnn/GSM8332741/suppl/GSM8332741_K562_2_chimericReads.csv.gz'

def _read_prim_seq_chimeric_reads(link: str, **kwargs) -> pd.DataFrame:
    header = pd.read_csv(
        link,
        nrows=1,
        header=None
    )
    result = _read_tsv(
        link,
        header=None,
        skiprows=1,
        names=list(header.iloc[0]) + ['R1GeneType', 'R2GeneType'],
        sep=',',
        **kwargs
    )
    return result


def _get_prim_seq_ids_mapping() -> pd.Series:
    data1 = _read_prim_seq_chimeric_reads(
        LINK1,
        usecols=['R1Tx', 'R1Gene', 'R2Tx', 'R2Gene']
    )
    data2 = _read_prim_seq_chimeric_reads(
        LINK2,
        usecols=['R1Tx', 'R1Gene', 'R2Tx', 'R2Gene']
    )
    data = pd.concat([
        data1[['R1Tx', 'R1Gene']].rename(columns={'R1Tx': 'Tx', 'R1Gene': 'Gene'}),
        data1[['R2Tx', 'R2Gene']].rename(columns={'R2Tx': 'Tx', 'R2Gene': 'Gene'}),
        data2[['R1Tx', 'R1Gene']].rename(columns={'R1Tx': 'Tx', 'R1Gene': 'Gene'}),
        data2[['R2Tx', 'R2Gene']].rename(columns={'R2Tx': 'Tx', 'R2Gene': 'Gene'}),
    ]).drop_duplicates()

    data['Tx'] = data['Tx'].str.split('.').str[0]
    data = data[data['Tx'].str[:2].isin(('NM', 'NR'))].copy()

    data['yagid'] = id2yagid(data['Tx'])

    result = data.groupby('Gene')['yagid'].min()

    return result

def load_prim_seq_data():
    response = requests.get('https://sysbiocomp.ucsd.edu/prim/HuRPA.csv', verify=False)
    response.raise_for_status()
    result = pd.read_csv(StringIO(response.text))

    gene2yagid = _get_prim_seq_ids_mapping()

    result['yagid'] = result['RNA'].map(gene2yagid)

    result['yapid'] = id2yapid('SYMBOL:' + result['protein'])

    result = result[result['yapid'].str.startswith('YAPID') & result['yagid'].str.startswith('YAGID')]

    result = result[['yagid', 'yapid']].drop_duplicates()

    return result
