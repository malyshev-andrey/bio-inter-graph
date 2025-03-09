import pandas as pd

from ..shared import memory, BED_COLUMNS
from .main import sanitize_bed


def _collapse_SPIN_states(states: pd.Series) -> pd.Series:
    result = states.replace(
        '^Near_Lm.*', 'Near_Lm', regex=True
    ).replace(
        '^Interior_Act.*', 'Interior_Act', regex=True
    ).replace(
        '^Interior_Repr.*', 'Interior_Repr', regex=True
    )
    return result


def _load_SPIN_annotation(**kwargs) -> pd.DataFrame:
    id = '1gdwtrhTctddO9TCBXBaZpZFOAHWCUTli'
    url = f'https://drive.usercontent.google.com/download?id={id}&export=download&confirm=t'

    default_kwargs = dict(
        sep='\t',
        header=None,
        dtype='str',
        names=BED_COLUMNS[:4]
    )
    default_kwargs.update(kwargs)

    result = pd.read_csv(url, **default_kwargs)

    result['start'] = result['start'].astype('int')
    result['end'] = result['end'].astype('int')
    assert (result['start'] < result['end']).all()

    result = result.rename(columns={'name': 'SPIN_full'})
    result['SPIN'] = _collapse_SPIN_states(result['SPIN_full'])

    return result


def _collapse_ChromHMM_states(states: pd.Series) -> pd.Series:
    result = states.replace(
        r'^Enh.*', 'Enh', regex=True
    ).replace(
        r'^Tss.*', 'Tss', regex=True
    ).replace(
        r'^Tx.*', 'Tx', regex=True
    ).replace(
        r'^ReprPC.*', 'ReprPC', regex=True
    ).replace(
        'ZNF/Rpts', 'Quies'
    )
    expected_states = ['Enh', 'ReprPC', 'Tss', 'Tx', 'Quies', 'Het']
    assert result.isin(expected_states).all()
    return result


@memory.cache
def load_ChromHMM_annotation(**kwargs):
    url = 'https://personal.broadinstitute.org/cboix/epimap/ChromHMM/observed_aux_18_hg38/CALLS/BSS00762_18_CALLS_segments.bed.gz'
    default_kwargs = dict(
        sep='\t',
        header=None,
        usecols=range(6),
        names=BED_COLUMNS,
        dtype='str'
    )
    default_kwargs.update(kwargs)
    result = pd.read_csv(url, **default_kwargs)
    result = sanitize_bed(result, stranded=False)
    assert (result['start'] < result['end']).all()

    assert result['score'].eq('0').all()
    assert result['strand'].eq('.').all()
    result = result.drop(columns=['score', 'strand'])

    columns_map = {f'{c}1': c for c in result.columns}

    n = result.shape[0]

    SPIN = _load_SPIN_annotation()
    result = ba.bed_intersect(
        result, SPIN,
        strandedness=None,
        unify_chr_assembly='hg38',
        jaccard=True,
        how='left'
    )
    result = result.rename(columns=columns_map)

    result = result.sort_values('jaccard', ascending=False)
    result = result.drop_duplicates(columns_map.values(), keep='first')
    assert result.shape[0] == n

    result = result.drop(columns=['start2', 'end2', 'Overlap', 'jaccard'])

    result = result.rename(columns={'name': 'state_full'})
    result['name'] = 'YALID' + result.index.astype('str').str.zfill(7)
    assert result['name'].str.len().eq(12).all()

    result['state'] = _collapse_ChromHMM_states(result['state_full'])

    return result
