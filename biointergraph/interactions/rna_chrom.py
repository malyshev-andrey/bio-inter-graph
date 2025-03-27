import pandas as pd

from ..shared import GOOGLE_DRIVE_URL, memory, _read_tsv
from ..annotations import load_extended_annotation, load_chromhmm_annotation
from ..ids_mapping import id2yagid
from .main import _annotate_peaks


@memory.cache
def load_redc_redchip_data() -> pd.DataFrame:
    result = _read_tsv(
        GOOGLE_DRIVE_URL.format(id='1nkg0Iofz8azz6BEfISWXG_DNQMlHEbi6'),
        compression='zip',
        desc='Red-C, RedChIP',
        usecols=lambda name: name not in {'ID', 'sample'}
    )

    names_map = {
        'seqid': 'chr1', 'start': 'gene_start',
        'end': 'gene_end', 'strand': 'strand1',
        'extended_gene_id': 'name'
    }
    annotation = load_extended_annotation()
    annotation = annotation.rename(columns=names_map)

    result = result.merge(annotation, how='left', validate='many_to_one')
    assert not result['name'].isna().any()

    cols = ['chr', 'start', 'end', 'name']
    result = result.rename(columns={f'{c}2': c for c in cols})
    result = result[cols]

    result = _annotate_peaks(
        result, load_chromhmm_annotation(),
        assembly='hg38', stranded=False,
        desc='Red-C, RedChIP'
    )

    result['source'] = id2yagid(result['source'], strict=True)
    result = result.drop_duplicates()

    return result
