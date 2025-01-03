from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
from tqdm.auto import tqdm

from ..shared import memory
from ..interactions.karr_seq import _retrieve_karr_seq_metadata, _load_single_karr_seq
from ..ids import drop_id_version


def _load_karr_seq_ids(path: str) -> pd.DataFrame:
    result = _load_single_karr_seq(
        path,
        filter_func=lambda df: df[['seqid1', 'seqid2']]
    )
    result = set(result['seqid1']) | set(result['seqid2'])
    return result


@memory.cache
def karr_seq_ids2entrezgene_id():
    ids = set()

    with ProcessPoolExecutor(max_workers=5) as executor:
        futures = []
        for url in _retrieve_karr_seq_metadata()['url']:
            futures.append(executor.submit(_load_karr_seq_ids, url,))

        for future in tqdm(as_completed(futures)):
            ids.update(future.result())


    ids = pd.Series(list(ids))
    ids = drop_id_version(ids)
    assert ids.is_unique
    return ids
