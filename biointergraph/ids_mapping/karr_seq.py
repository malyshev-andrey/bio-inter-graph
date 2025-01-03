import pandas as pd

from ..shared import memory
from ..interactions.karr_seq import _retrieve_karr_seq_metadata, _load_single_karr_seq
from ..ids import drop_id_version


@memory.cache
def karr_seq_ids2entrezgene_id():
    ids = set()
    for url in _retrieve_karr_seq_metadata()['url']:
        data = _load_single_karr_seq(url, filter_func=lambda df: df[['seqid1', 'seqid2']])
        ids.update(drop_id_version(data['seqid1']))
        ids.update(drop_id_version(data['seqid2']))

    return pd.Series(ids)
