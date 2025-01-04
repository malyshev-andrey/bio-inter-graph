from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from tqdm.auto import tqdm

from .entrez import refseq_transcript_id2entrez_gene_id
from ..shared import memory
from ..interactions.karr_seq import _retrieve_karr_seq_metadata, _load_single_karr_seq
from ..ids import drop_id_version


@memory.cache
def karr_seq_ids2entrezgene_id():
    ids = set()

    with ThreadPoolExecutor(max_workers=20) as executor:
        futures = []
        for url in _retrieve_karr_seq_metadata()['url']:
            futures.append(executor.submit(
                _load_single_karr_seq,
                url,
                filter_func=lambda df: df[['seqid1', 'seqid2']]
            ))

        for future in as_completed(futures):
            result = future.result()
            ids.update(result['seqid1'])
            ids.update(result['seqid2'])


    ids = pd.Series(list(ids))
    ids = drop_id_version(ids)
    assert ids.is_unique
    result = pd.concat(
        [ids, refseq_transcript_id2entrez_gene_id(ids)],
        axis=1
    )
    result.columns = 'refseq_transcript_id', 'entrezgene_id'
    return result
