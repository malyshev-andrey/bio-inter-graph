import pandas as pd

from ..shared import memory, BED_COLUMNS
from .main import sanitize_bed


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
def load_ChromHMM_annotation() -> pd.DataFrame:
    url = 'https://personal.broadinstitute.org/cboix/epimap/ChromHMM/observed_aux_18_hg38/CALLS/BSS00762_18_CALLS_segments.bed.gz'
    result = pd.read_csv(url, sep='\t', header=None, usecols=range(6), names=BED_COLUMNS, dtype='str')
    result = sanitize_bed(result, stranded=False)
    assert (result['start'] < result['end']).all()

    result = result.rename(columns={'name': 'state_full'})
    result['name'] = 'YALID' + result.index.astype('str').str.zfill(7)
    assert result['name'].str.len().eq(12).all()

    result['state'] = _collapse_ChromHMM_states(result['state_full'])
    return result
