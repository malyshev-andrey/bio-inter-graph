import importlib.resources

import pandas as pd

from ..shared import memory, BED_COLUMNS, GOOGLE_DRIVE_URL, _read_tsv, _df_hash
from .main import sanitize_bed, _split_annotation_into_bins
from .intersect import best_left_intersect


def _merge_spin_states(states: pd.Series) -> pd.Series:
    states = states.case_when([
        (states.str.match('^Near_Lm.*$'), 'Near_Lm'),
        (states.str.match('^Interior_Act.*$'), 'Interior_Act'),
        (states.str.match('^Interior_Repr.*$'), 'Interior_Repr')
    ])

    expected_states = [
        'Interior_Act', 'Interior_Repr', 'Lamina',
        'Lamina_Like', 'Near_Lm', 'Speckle'
    ]
    assert states.isin(expected_states).all()
    assert states.nunique(dropna=False) == len(expected_states)
    assert not states.isna().any()

    return states


def _load_spin_annotation() -> pd.DataFrame:
    result = _read_tsv(
        GOOGLE_DRIVE_URL.format(id='1gdwtrhTctddO9TCBXBaZpZFOAHWCUTli'),
        header=None,
        names=BED_COLUMNS[:4],
        chunksize=None
    )
    result = sanitize_bed(result)

    result = result.rename(columns={'name': 'SPIN_full'})
    result['SPIN'] = _merge_spin_states(result['SPIN_full'])

    return result


def _merge_chromhmm_states(states: pd.Series) -> pd.Series:
    states = states.case_when([
        (states.str.match(r'^Enh.*$'), 'Enh'),
        (states.str.match(r'^Tss.*$'), 'Tss'),
        (states.str.match(r'^Tx.*$'), 'Tx'),
        (states.str.match(r'^ReprPC.*$'), 'ReprPC'),
        (states.eq('ZNF/Rpts'), 'Quies')
    ])

    expected_states = ['Enh', 'ReprPC', 'Tss', 'Tx', 'Quies', 'Het']
    assert states.isin(expected_states).all()
    assert states.nunique(dropna=False) == len(expected_states)
    assert not states.isna().any()

    return states


@memory.cache
def load_chromhmm_annotation(split_bin: int|None = None) -> pd.DataFrame:
    REBUILD_CHROMHMM_ANNOTATION = False

    if split_bin == 500 and not REBUILD_CHROMHMM_ANNOTATION:
        with importlib.resources.open_binary('bio-inter-graph.static', 'chromhmm_500.tsv.gz') as file:
            result = pd.read_csv(
                file, compression='gzip',
                sep='\t', header=None
            )

        assert _df_hash(result) == 'e1f53aa30ad4c3303ae55fc0b5430daf8b8e379f'

        return result

    result = _read_tsv(
        'https://personal.broadinstitute.org/cboix/epimap/ChromHMM/observed_aux_18_hg38/CALLS/BSS00762_18_CALLS_segments.bed.gz',
        header=None,
        usecols=range(6),
        names=BED_COLUMNS
    )
    result = sanitize_bed(result, stranded=False)
    result = result.drop(columns=['score', 'strand'])

    result = result.rename(columns={'name': 'state_full'})
    result['state'] = _merge_chromhmm_states(result['state_full'])

    if split_bin is not None:
        result = _split_annotation_into_bins(result, bin_size=split_bin)

    result = best_left_intersect(
        result,
        _load_spin_annotation(),
        stranded=False,
        unify_chr_assembly='hg38'
    )

    result = result.drop(columns=['start2', 'end2', 'jaccard'])

    result = result.sort_values(['chr', 'start', 'end'])
    result = result.reset_index(drop=True)
    result['name'] = 'YALID' + result.index.astype('str').str.zfill(7)
    assert result['name'].str.len().eq(12).all()

    if split_bin == 500:
        assert _df_hash(result) == 'e1f53aa30ad4c3303ae55fc0b5430daf8b8e379f'

    return result


def yalid2state(ids: pd.Series|None = None) -> pd.Series:
    result = load_chromhmm_annotation()
    result = result.set_index('name')
    result = result['state']

    if ids is not None:
        result = ids.map(result)

    return result
