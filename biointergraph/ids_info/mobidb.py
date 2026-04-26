import requests

import pandas as pd
from tqdm.auto import tqdm

from ..shared import memory
from ..ids_mapping import yapid2ids_by_type
from .uniprot import uniprot_id_info


MOBIDB_API_URL = 'https://mobidb.org/api/download_page'
BATCH_SIZE = 200

PREDICTED_FEATURE = 'prediction-disorder-mobidb_lite'
CURATED_FEATURE = 'curated-disorder-disprot'


@memory.cache
def mobidb_disorder_info(organism_id: str = '9606') -> pd.DataFrame:
    human_accessions = uniprot_id_info(organism_id=organism_id)
    human_accessions = set(human_accessions['Entry'])

    batches = sorted(human_accessions)
    batches = [
        batches[i:i + BATCH_SIZE]
        for i in range(0, len(batches), BATCH_SIZE)
    ]

    records = []
    for batch in tqdm(batches, desc='MobiDB disorder'):
        response = requests.get(MOBIDB_API_URL, params={
            'acc': ','.join(batch),
            'format': 'json'
        })
        response.raise_for_status()
        data = response.json()['data']

        for entry in data:
            acc = entry.get('acc')
            organism = entry.get('organism', '')

            if 'Homo sapiens' not in organism:
                continue

            predicted = entry.get(PREDICTED_FEATURE, {})
            curated = entry.get(CURATED_FEATURE, {})

            records.append({
                'acc': acc,
                'disorder_predicted': predicted.get('content_fraction', 0.0),
                'disorder_curated': curated.get('content_fraction', float('nan')),
            })

    result = pd.DataFrame(records)

    if result.empty:
        result = pd.DataFrame(columns=['acc', 'disorder_predicted', 'disorder_curated'])

    result = result.set_index('acc', verify_integrity=True)
    result['disorder_predicted'] = result['disorder_predicted'].astype(float)
    result['disorder_curated'] = result['disorder_curated'].astype(float)

    assert result['disorder_predicted'].between(0, 1).all()
    assert result['disorder_curated'].dropna().between(0, 1).all()

    return result


def yapid2is_disordered(curated_only: bool = False) -> pd.Series:
    disorder_info = mobidb_disorder_info()

    col = 'disorder_curated' if curated_only else 'disorder_predicted'
    uniprot2disorder = disorder_info[col]

    result = yapid2ids_by_type()['uniprot'].explode()
    result = result.map(uniprot2disorder).fillna(0)
    result = result.groupby(level=0).max()

    return result
