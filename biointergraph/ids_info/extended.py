import pandas as pd

from biointergraph.shared import UNIFY_BIOTYPES
from ..annotations import load_extended_annotation


def extended_gene_id2biotype(ids: pd.Series|None = None) -> pd.Series:
    result = load_extended_annotation()
    result = result.set_index('extended_gene_id', verify_integrity=True)
    result = result['gene_type']

    result = result.where(
        ~result.str.contains('lncRNA'),
        'lncRNA'
    )
    result = result.where(
        ~result.str.contains('pseudogene'),
        'pseudogene'
    )
    result = result.replace({
        'uncertain_coding': 'mRNA',
        'vlinc': 'lncRNA',
        'small_RNA': 'sRNA',
        'trna': 'tRNA',
        'scaRna': 'scaRNA',
        'short_ncRNA': float('nan'),
        'RNA': float('nan')
    })
    result = result.replace(UNIFY_BIOTYPES)

    if ids is not None:
        result = ids.map(result)
    return result
