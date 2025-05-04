from itertools import combinations
from concurrent.futures import ThreadPoolExecutor, as_completed
from importlib.resources import files

import pandas as pd
import networkx as nx

from .intersect import (
    gencode_refseq_intersect2pairwise,
    extended_refseq_intersect2pairwise,
    extended_gencode_intersect2pairwise
)
from .BioMart import load_BioMart_pairwise
from .OrgHsEgDb import load_OrgHsEgDb_pairwise
from .entrez import karr_seq_ids2entrezgene_id
from ..annotations import extended_gene_id2ensembl_gene_id, load_extended_annotation
from ..shared import ID_TYPES, memory
from ..ids import drop_id_version


@memory.cache
def _build_yagid_graph() -> pd.Series:
    REBUILD_YAGID_MAPPING = False

    if not REBUILD_YAGID_MAPPING:
        with (files('biointergraph.static') / "id2yagid.json").open('r') as file:
            result = pd.read_json(file, typ='series')
        return result


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
    pairs.append(extended_gencode_intersect2pairwise())
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

    result = {}
    for i, component in enumerate(nx.connected_components(yagid_graph)):
        assert i < 1e7
        for node in component:
            result[node] = 'YAGID' + str(i).zfill(7)
    result = pd.Series(result)

    n_components = result.nunique()
    print(f'Build YAGID graph: {len(result)} ids, {n_components} components')
    is_valid = result.index.str.match(r'^\d+|ENS[TG]\d{11}|N[MR]_\d+|EXTG\d{7}$')
    sample = result[~is_valid].index[:5]
    assert is_valid.all(), f'_build_yagid_graph: invalid IDs: {sample}'
    return result


def id2yagid(ids: pd.Series|None = None, *, strict: bool = False) -> pd.Series:
    result = _build_yagid_graph()
    result.name = 'yagid'

    if ids is not None:
        ids = drop_id_version(ids)
        result = ids.map(result)
        if strict:
            assert not result.isna().any()
        result = result.combine_first(ids)
    return result


def yagid2ids(yagid: str|list[str]|None = None, *, squeeze: bool = True) -> pd.Series|list[str]:
    result = id2yagid()
    if yagid is not None:
        if isinstance(yagid, list):
            result = result[result.isin(yagid)]
        else:
            assert isinstance(yagid, str)
            result = result[result.eq(yagid)]
    result = result.to_frame().reset_index(names='ids')
    result = result.groupby('yagid')['ids'].agg(list)

    if squeeze and isinstance(yagid, str):
        result = result.item()
    return result
