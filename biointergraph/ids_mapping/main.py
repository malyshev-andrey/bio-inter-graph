from itertools import combinations
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import networkx as nx

from .BioMart import load_BioMart_pairwise
from .OrgHsEgDb import load_OrgHsEgDb_pairwise
from ..shared import ID_TYPES, memory
from ..annotations import gencode_refseq_intersect2pairwise


@memory.cache
def _build_yagid_graph():
    pairs = []
    # BioMart + OrgHsEgDb data
    with ThreadPoolExecutor(max_workers=6) as executor:
        futures = []
        for ids in combinations(ID_TYPES, r=2):
            for assembly in 'GRCh37', 'GRCh38':
                futures.append(executor.submit(load_BioMart_pairwise, *ids, assembly=assembly))

            futures.append(executor.submit(load_OrgHsEgDb_pairwise, *ids))

        for future in as_completed(futures):
            pairs.append(future.result())

    # annotations intersections
    pairs.append(gencode_refseq_intersect2pairwise('hg19'))
    pairs.append(gencode_refseq_intersect2pairwise('hg38'))

    for df in pairs:
        assert df.shape[1] == 2
        df.columns = 'source', 'target'
    pairs = pd.concat(pairs)
    assert pairs.shape[1] == 2

    yagid_graph = nx.from_pandas_edgelist(pairs)
    connected_components = nx.connected_components(yagid_graph)
    result = {}
    for i, component in enumerate(connected_components):
        assert i < 1e6
        for node in component:
            result[node] = 'YAGID' + str(i).zfill(6)
    result = pd.Series(result)

    n_components = result.nunique()
    print(f'Build YAGID graph: {len(result)} ids, {n_components} components')
    is_valid = result.index.str.match(r'^\d+|ENS[TG]\d{11}|N[MR]_\d+$')
    sample = result[~is_valid].index[:5]
    assert is_valid.all(), f'_build_yagid_graph: invalid IDs: {sample}'
    return result


def id2yagid(ids: pd.Series) -> pd.Series:
    mapping = _build_yagid_graph()
    ids = ids.map(mapping)
    return ids
