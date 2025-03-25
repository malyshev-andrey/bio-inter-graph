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
    result = result.where(
        ~result.str.match('^(TR|IG)_[A-Z]_gene$'),
        'IG_TR_gene'
    )
    result = result.replace({
        'uncertain_coding': 'mRNA',
        'vlinc': 'lncRNA',
        'small_RNA': 'sRNA',
        'trna': 'tRNA',
        'scaRna': 'scaRNA',
        'HAcaBox': 'snoRNA',
        'CDBox': 'snoRNA'
    })
    result = result.replace(
        ['short_ncRNA', 'RNA', 'sense_overlap_RNA', 'structural_RNA', '__na'],
        float('nan')
    )
    result = result.replace(UNIFY_BIOTYPES)

    if ids is not None:
        result = ids.map(result)
    return result
