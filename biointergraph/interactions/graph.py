from concurrent.futures import ThreadPoolExecutor, as_completed
from time import sleep
import importlib.resources
from typing import Callable

import requests
import pandas as pd
import networkx as nx
from tqdm.auto import tqdm

from biointergraph.ids_mapping.main import id2yagid

from .encode import (
    load_encode_eclip_data,
    load_encode_rip_data,
    load_encode_iclip_data,
    load_encode_chip_seq_data
)
from .rna_protein import load_postar3_data, load_frip_seq_data
from .karr_seq import load_karr_seq_data
from .ric_seq import load_ric_seq_data
from .protein import load_intact_interactions, load_biogrid_interactions, load_string_interactions
from .rna_chrom import load_redc_redchip_data
from .gtrd import load_gtrd_chip_seq_data
from ..shared import memory
from ..ids_mapping.protein import _build_yapid_graph, id2yapid
from ..annotations import yalid2state
from ..ids_mapping import id2yagid, id2yapid
from ..ids_info import yagid2biotype


def _wrapper(dataset: str, func: Callable, **kwargs) -> pd.DataFrame:
    result = func(**kwargs)
    result = result.reset_index(drop=True)
    assert result.shape[1] == 2
    result.columns = 'source', 'target'
    assert (result['source'] != result['target']).all()

    n_pairs = result.shape[0]
    swap_mask = result['source'] > result['target']
    result.loc[swap_mask, ['source', 'target']] = result.loc[swap_mask, ['target', 'source']].values
    assert (result['source'] < result['target']).all()
    assert not result.duplicated().any()

    result['dataset'] = dataset
    assert result.shape[0] == n_pairs

    return result


def _remove_minor_components(graph: nx.Graph) -> nx.Graph:
    sizes = [len(c) for c in nx.connected_components(graph)]
    if len(sizes) < 2:
        return graph

    sizes.sort()
    print(f'The largest graph components: {sizes[-1]}, {sizes[-2]}, {sizes[-3]}')
    major_component_size = sizes[-1]

    nodes_to_delete = []
    for component in nx.connected_components(graph):
        if len(component) < major_component_size:
            nodes_to_delete.extend(component)
    graph.remove_nodes_from(nodes_to_delete)
    assert nx.number_connected_components(graph) == 1

    print(f'Minor components nodes removed: {len(nodes_to_delete)}')
    return graph


@memory.cache
def build_main_graph(max_workers: int = 2) -> nx.Graph:
    REBUILD_MAIN_GRAPH = False

    if not REBUILD_MAIN_GRAPH:
        with importlib.resources.open_binary('bio-inter-graph.static', 'edges.tsv.gz') as file:
            result = pd.read_csv(file, compression='gzip', sep='\t', dtype='str')

        result = nx.from_pandas_edgelist(result, edge_attr='dataset')

        return result


    data = [
        ('KARR-seq', load_karr_seq_data, dict(cell_line='K562', pvalue=0.05)),
        ('ENCODE eCLIP', load_encode_eclip_data, dict(assembly='hg38', annotation='gencode', cell_line='K562')),
        ('ENCODE RIP', load_encode_rip_data, dict(annotation='gencode', cell_line='K562')),
        ('ENCODE iCLIP', load_encode_iclip_data, dict(annotation='gencode', cell_line='K562')),
        ('POSTAR3', load_postar3_data, dict(species='human', cell_line='K562', annotation='gencode')),
        ('fRIP-seq', load_frip_seq_data, dict()),
        ('RIC-seq', load_ric_seq_data, dict(pvalue=0.05)),
        ('IntAct', load_intact_interactions, dict()),
        ('BioGRID', load_biogrid_interactions, dict()),
        ('STRING', load_string_interactions, dict(min_score=700)),
        ('ENCODE ChIP-seq', load_encode_chip_seq_data, dict(assembly='hg38', cell_line='K562')),
        ('Red-C & RedChIP', load_redc_redchip_data, dict()),
        ('GTRD', load_gtrd_chip_seq_data, dict(cell_line='K562'))
    ]

    tqdm_kwargs = dict(total=len(data), unit='dataset', desc='Collecting data: ')
    if max_workers > 1:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = []
            for dataset, func, kwargs in data:
                futures.append(executor.submit(_wrapper, dataset, func, **kwargs))
                sleep(1)

            data = [f.result() for f in tqdm(as_completed(futures), **tqdm_kwargs)]

    else:
        data = [_wrapper(dataset, func, **kwargs) for dataset, func, kwargs in tqdm(data, **tqdm_kwargs)]

    data = pd.concat(data)
    assert not data.duplicated().any()

    assert data.shape[1] == 3
    regex = r'^YA[LPG]ID\d{7}$'
    assert data['source'].str.match(regex).all()
    assert data['target'].str.match(regex).all()

    data = data.groupby(['source', 'target'], as_index=False, observed=True)['dataset'].agg(','.join)

    graph = nx.from_pandas_edgelist(data, edge_attr='dataset')

    graph = _remove_minor_components(graph)

    return graph


def _node_id2node_type(ids: pd.Series) -> pd.Series:
    prefix_map = {
        'YAPID': 'protein',
        'YAGID': 'RNA',
        'YALID': 'DNA'
    }
    result = ids.str[:5].map(prefix_map)
    assert not result.isna().any()

    return result


def _describe_nodes(graph: nx.Graph) -> pd.DataFrame:
    result = pd.DataFrame(graph.degree(), columns=['node', 'degree'])
    assert result['node'].str.match(r'^YA[GPL]ID\d{7}$').all()

    result['type'] = _node_id2node_type(result['node'])

    return result


def _symmetric_crosstab(pairs: pd.DataFrame) -> pd.DataFrame:
    assert pairs.shape[1] == 2

    c1, c2 = pairs.columns
    pairs = pd.concat([
        pairs,
        pairs[pairs[c1] != pairs[c2]].rename(columns={c1: c2, c2: c1})
    ])

    result = pd.crosstab(pairs[c1], pairs[c2])
    return result


def describe_graph(graph: nx.Graph) -> None:
    nodes = _describe_nodes(graph)
    print(f'Nodes: {nodes.shape[0]}')
    print(pd.DataFrame({
        'count': nodes['type'].value_counts(),
        'frac':nodes['type'].value_counts(normalize=True)
    }))

    edges = pd.DataFrame(graph.edges(), columns=['source', 'target'])
    print(f'Edges: {edges.shape[0]}')
    edges['source_type'] = _node_id2node_type(edges['source'])
    edges['target_type'] = _node_id2node_type(edges['target'])
    print(_symmetric_crosstab(edges[['source_type', 'target_type']]))

    print('Degrees distribution:')
    stats = nodes['degree'].describe()
    print('\t' + str(stats[:10]).replace('\n', '\n\t'))

    iterator = tqdm(range(10)) if nodes.shape[0] > 10000 else range(10)
    diameter_sample = [nx.approximation.diameter(graph) for _ in iterator]
    print(f'Diameter: {min(diameter_sample)}-{max(diameter_sample)}')


def _community2enrichment(nodes: list[str]) -> tuple[str, float]:
    nodes = [n for n in nodes if n.startswith('YAPID')]
    if not nodes:
        return '', 1

    id2yapid = _build_yapid_graph()
    ids = pd.Series(id2yapid[id2yapid.isin(nodes)].index)
    ids = ids.str.removeprefix('SYMBOL:')

    response = requests.post(
        url = 'https://biit.cs.ut.ee/gprofiler/api/gost/profile/',
        json = {
            'organism': 'hsapiens',
            'query': ids.tolist(),
            'sources': ['GO'],
            'numeric_ns': 'BIOGRID'
        }
    )
    result = pd.DataFrame(response.json()['result'])

    if result.shape[0] == 0:
        return '', 1

    result = result.iloc[0]
    return result['name'], result['p_value']


def detect_communities(graph) -> pd.DataFrame:
    communities = nx.community.louvain_communities(graph, resolution=2, seed=42)
    result = []
    for community in communities:
        community = list(community)
        result.append({
            'nodes': community,
            'size': len(community),
            'protein_frac': pd.Series(community).str.startswith('YAPID').mean()
        })
    result = pd.DataFrame(result)

    tqdm.pandas()
    go_terms = result['nodes'].progress_apply(_community2enrichment)
    go_terms = go_terms.tolist()
    result[['go_term', 'p_value']] = pd.DataFrame(go_terms)

    result = result.sort_values('p_value')

    return result


def _node2neighbors_types(graph, binary: bool = False) -> pd.DataFrame:
    edges = pd.DataFrame(graph.edges(), columns=['source', 'target'])

    swap = {'target': 'source', 'source': 'target'}
    edges = pd.concat([edges, edges.rename(columns=swap)])

    self_loops = edges.drop_duplicates('source').copy()
    self_loops['target'] = self_loops['source']

    edges = pd.concat([edges, self_loops])

    edges['target'] = _node_id2node_type(edges['target'])

    result = pd.crosstab(edges['source'], edges['target'])
    if binary:
        result = result > 0

    result = result.reset_index(names='node')

    return result


def describe_nodes(graph: nx.Graph) -> pd.DataFrame:
    result = _describe_nodes(graph)

    ids_mapping = pd.concat([id2yagid(), id2yapid()])
    ids_mapping.name = 'node'
    ids_mapping = ids_mapping.to_frame().reset_index(names='id')

    ids_mapping = ids_mapping.groupby('node')['id'].agg(pd.Series.to_list)

    result['ids'] = result['node'].map(ids_mapping)

    tqdm.pandas()
    result['neighbors'] = result['node'].progress_apply(graph.neighbors)
    result['neighbors'] = result['neighbors'].progress_apply(list)

    result['subtype'] = result['type'].case_when([
        (result['type'].eq('RNA'), yagid2biotype(result['node'])),
        (result['type'].eq('DNA'), yalid2state(result['node']))
    ])

    result = result.merge(
        _node2neighbors_types(graph),
        how='left',
        on='node',
        validate='one_to_one'
    )

    assert (result[['DNA', 'protein', 'RNA']].sum(axis=1) == (result['degree'] + 1)).all()

    return result


def id2subgraph(graph: nx.Graph, id: str) -> nx.Graph:
    yagid = id2yagid()[id]
    neighbors = list(graph.neighbors(yagid))
    result = graph.subgraph(neighbors + [yagid])
    return result


def _lighten_graph(graph: nx.Graph, *, inplace: bool = False) -> nx.Graph:
    if not inplace:
        graph = graph.copy()

    nodes = describe_nodes(graph)
    nodes_to_remove = (
        (nodes['degree'].eq(1) & nodes['type'].eq('DNA'))
        | nodes['subtype'].eq('mRNA')
    )
    nodes_to_remove = nodes.loc[nodes_to_remove, 'node']
    graph.remove_nodes_from(nodes_to_remove)
    return graph


def build_light_graph(max_workers: int = 2) -> nx.Graph:
    graph = build_main_graph(max_workers=max_workers)
    _lighten_graph(graph, inplace=True)
    graph = _remove_minor_components(graph)
    return graph
