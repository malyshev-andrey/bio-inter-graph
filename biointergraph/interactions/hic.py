import pandas as pd

from ..shared import _read_tsv, GOOGLE_DRIVE_URL
from ..annotations import load_chromhmm_annotation, best_left_intersect

INTER_DATA_ID = '1nQ4c3nkcFCzr-rlnfxvItk5b_vStk4tq'
INTRA_DATA_ID = '11XvGPgC9FEF1VDoO-x6-nVCbnLpGc8q_'


def _get_bin2yalid_map(bins: pd.DataFrame) -> pd.DataFrame:
    result = best_left_intersect(
        bins,
        load_chromhmm_annotation(),
        stranded=False,
        unify_chr_assembly='hg38'
    )
    assert not bins[['chr', 'mid']].duplicated().any()
    return result[['chr', 'mid', 'name']]


def load_hic_data() -> pd.DataFrame:
    result = pd.concat([
        _read_tsv(
            GOOGLE_DRIVE_URL.format(id=INTER_DATA_ID),
            usecols=['chr1', 'fragmentMid1', 'chr2', 'fragmentMid2', 'q-value'],
            compression='gzip',
            dtype=None
        ),
        _read_tsv(
            GOOGLE_DRIVE_URL.format(id=INTRA_DATA_ID),
            compression='gzip',
            dtype=None
        )
    ])

    bins = pd.concat([
        result[['chr1', 'fragmentMid1']].rename(columns={'chr1': 'chr', 'fragmentMid1': 'mid'}),
        result[['chr2', 'fragmentMid2']].rename(columns={'chr2': 'chr', 'fragmentMid2': 'mid'})
    ]).drop_duplicates()

    bins = bins.assign(
        start=bins['mid']-500,
        end=bins['mid']+499
    )

    bin2yalid_map = _get_bin2yalid_map(bins)

    result = result.merge(
        bin2yalid_map.rename(columns={'chr': 'chr1', 'mid': 'fragmentMid1', 'name': 'yalid1'}),
        how='left',
        validate='many_to_one'
    ).merge(
        bin2yalid_map.rename(columns={'chr': 'chr2', 'mid': 'fragmentMid2', 'name': 'yalid2'}),
        how='left',
        validate='many_to_one'
    )

    result = result.groupby(['yalid1', 'yalid2'], as_index=False, observed=True)['q-value'].min()

    assert result.notna().all().all()

    return result
