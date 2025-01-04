from itertools import combinations
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import networkx as nx

from .intersect import gencode_refseq_intersect2pairwise, extended_refseq_intersect2pairwise
from .BioMart import load_BioMart_pairwise
from .OrgHsEgDb import load_OrgHsEgDb_pairwise
from .entrez import karr_seq_ids2entrezgene_id
from ..annotations import extended_gene_id2ensembl_gene_id, load_extended_annotation
from ..shared import ID_TYPES, memory


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

    # extended annotation
    pairs.append(extended_refseq_intersect2pairwise())
    pairs.append(extended_gene_id2ensembl_gene_id())

    # KARR-seq refseq_transcript_ids -> entrezgene_id
    pairs.append(karr_seq_ids2entrezgene_id())

    for df in pairs:
        assert df.shape[1] == 2
        df.columns = 'source', 'target'
    pairs = pd.concat(pairs)
    assert pairs.shape[1] == 2

    yagid_graph = nx.from_pandas_edgelist(pairs)
    yagid_graph.add_nodes_from(load_extended_annotation()['extended_gene_id'])

    connected_components = nx.connected_components(yagid_graph)
    result = {}
    for i, component in enumerate(connected_components):
        assert i < 1e6
        for node in component:
            result[node] = 'YAGID' + str(i).zfill(7)
    result = pd.Series(result)

    n_components = result.nunique()
    print(f'Build YAGID graph: {len(result)} ids, {n_components} components')
    is_valid = result.index.str.match(r'^\d+|ENS[TG]\d{11}|N[MR]_\d+|EXTG\d{7}$')
    sample = result[~is_valid].index[:5]
    assert is_valid.all(), f'_build_yagid_graph: invalid IDs: {sample}'
    return result


def id2yagid(ids: pd.Series) -> pd.Series:
    mapping = _build_yagid_graph()
    ids = ids.map(mapping).combine_first(ids)
    return ids
