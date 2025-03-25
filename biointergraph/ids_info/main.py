from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

from .ensembl import ensembl_transcript_id2biotype, ensembl_gene_id2biotype
from .entrez import entrezgene_id2biotype
from .refseq import refseq_transcript_id2biotype
from .extended import extended_gene_id2biotype
from ..ids_mapping import id2yagid


def yagid2biotype(ids: pd.Series|None = None, return_weights: bool = False) -> pd.Series:
    data = [
        ensembl_transcript_id2biotype,
        ensembl_gene_id2biotype,
        entrezgene_id2biotype,
        refseq_transcript_id2biotype,
        extended_gene_id2biotype
    ]

    with ThreadPoolExecutor(max_workers=len(data)) as executor:
        data = [executor.submit(f) for f in data]
        id2biotype = pd.concat(f.result() for f in as_completed(data))

    result = id2yagid().to_frame().reset_index(names='id')

    n = id2yagid().nunique()

    result['biotype'] = result['id'].map(id2biotype)

    prefix_map = {
        'ENST': 'ensembl_transcript_id',
        'EXTG': 'extended_gene_id',
        'ENSG': 'ensembl_gene_id'
    }
    result['id_type'] = result['id'].case_when([
        (result['id'].str.isdigit(), 'entrezgene_id'),
        (result['id'].str[:3].isin({'NM_', 'NR_'}), 'refseq_transcript_id'),
        (result['id'].str[:4].isin(prefix_map), result['id'].str[:4].map(prefix_map))
    ])
    assert result['id_type'].nunique(dropna=False) == 5

    result['weight'] = 1 / result.groupby(['id_type', 'yagid']).transform('size')
    result['weight'] = result['weight'].mask(result['biotype'].isna(), 0)

    result = result.groupby(['yagid', 'biotype'], as_index=False, dropna=False)['weight'].sum()

    result['weight'] /= result.groupby('yagid')['weight'].transform('sum')

    result = result.sort_values('weight')
    result = result.drop_duplicates('yagid', keep='last')

    assert result.shape[0] == n

    result['biotype'] = result['biotype'].where(result['weight'] > 0.5, float('nan'))

    result = result.set_index('yagid', verify_integrity=True)
    weights, result = result['weight'], result['biotype']

    if ids is not None:
        result = ids.map(result)
        weights = ids.map(weights)

    return (result, weights) if return_weights else result
