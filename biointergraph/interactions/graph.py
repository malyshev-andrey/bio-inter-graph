from concurrent.futures import ThreadPoolExecutor, as_completed
from time import sleep
from importlib.resources import files
from typing import Callable
import random

import requests
import pandas as pd
import networkx as nx
from tqdm.auto import tqdm
from scipy.stats import entropy

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
from .prim_seq import load_prim_seq_data
from ..shared import memory
from ..annotations import yalid2state
from ..ids_mapping import id2yagid, yagid2ids, yapid2ids, yapid2best_id
from ..ids_info import yagid2biotype, yapid2is_nuclear


def _wrapper(dataset: str, func: Callable, **kwargs) -> pd.DataFrame:
    result = func(**kwargs)
    result = result.reset_index(drop=True)
    assert result.shape[1] == 2
    result.columns = 'source', 'target'
    assert (result['source'] != result['target']).all()

    n_pairs = result.shape[0]
    swap_mask = result['source'] > result['target']
    result.loc[swap_mask, ['source', 'target']
               ] = result.loc[swap_mask, ['target', 'source']].values
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
    print(
        f'The largest graph components: {sizes[-1]}, {sizes[-2]}, {sizes[-3]}')
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
    REBUILD_MAIN_GRAPH = True

    if not REBUILD_MAIN_GRAPH:
        with (files('biointergraph.static') / "edges.tsv.gz").open('rb') as file:
            result = pd.read_csv(file, compression='gzip',
                                 sep='\t', dtype='str')

        result = nx.from_pandas_edgelist(result, edge_attr='dataset')

        return result

    data = [
        ('KARR-seq', load_karr_seq_data, dict(cell_line='K562', pvalue=0.05)),
        ('ENCODE eCLIP', load_encode_eclip_data, dict(
            assembly='hg38', annotation='gencode', cell_line='K562')),
        ('ENCODE RIP', load_encode_rip_data, dict(
            annotation='gencode', cell_line='K562')),
        ('ENCODE iCLIP', load_encode_iclip_data, dict(
            annotation='gencode', cell_line='K562')),
        ('POSTAR3', load_postar3_data, dict(
            species='human', cell_line='K562', annotation='gencode')),
        ('fRIP-seq', load_frip_seq_data, dict()),
        ('RIC-seq', load_ric_seq_data, dict(pvalue=0.05)),
        ('IntAct', load_intact_interactions, dict()),
        ('BioGRID', load_biogrid_interactions, dict()),
        ('STRING', load_string_interactions, dict(min_score=700)),
        ('ENCODE ChIP-seq', load_encode_chip_seq_data,
         dict(assembly='hg38', cell_line='K562')),
        ('Red-C & RedChIP', load_redc_redchip_data, dict()),
        ('GTRD', load_gtrd_chip_seq_data, dict(cell_line='K562')),
        ('PRIM-seq', load_prim_seq_data, dict())
    ]

    tqdm_kwargs = dict(total=len(data), unit='dataset',
                       desc='Collecting data: ')
    if max_workers > 1:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = []
            for dataset, func, kwargs in data:
                futures.append(executor.submit(
                    _wrapper, dataset, func, **kwargs))
                sleep(1)

            data = [f.result()
                    for f in tqdm(as_completed(futures), **tqdm_kwargs)]

    else:
        data = [_wrapper(dataset, func, **kwargs)
                for dataset, func, kwargs in tqdm(data, **tqdm_kwargs)]

    data = pd.concat(data)
    assert not data.duplicated().any()

    assert data.shape[1] == 3
    regex = r'^YA[LPG]ID\d{7}$'
    assert data['source'].str.match(regex).all()
    assert data['target'].str.match(regex).all()

    data = data.groupby(['source', 'target'], as_index=False, observed=True)[
        'dataset'].agg(','.join)

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
        'frac': nodes['type'].value_counts(normalize=True)
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


def describe_nodes(
    graph: nx.Graph, *,
    subtypes: bool = True,
    neighbors_types: bool = True
) -> pd.DataFrame:
    result = _describe_nodes(graph)

    result['ids'] = result['node'].map(pd.concat([yagid2ids(), yapid2ids()]))

    tqdm.pandas(desc='Nodes processing: ')
    result['neighbors'] = result['node'].progress_apply(graph.neighbors)
    result['neighbors'] = result['neighbors'].progress_apply(list)

    if subtypes:
        result['subtype'] = result['type'].case_when([
            (result['type'].eq('RNA'), yagid2biotype(result['node'])),
            (result['type'].eq('DNA'), yalid2state(result['node']))
        ])

    if neighbors_types:
        result = result.merge(
            _node2neighbors_types(graph),
            how='left',
            on='node',
            validate='one_to_one'
        )
        assert (result[['DNA', 'protein', 'RNA']].sum(
            axis=1) - result['degree']).eq(1).all()

    return result


def id2subgraph(graph: nx.Graph, id: str) -> nx.Graph:
    yagid = id2yagid()[id]
    neighbors = list(graph.neighbors(yagid))
    result = graph.subgraph(neighbors + [yagid])
    return result


def _lighten_graph(graph: nx.Graph, *, inplace: bool = False) -> nx.Graph:
    if not inplace:
        graph = graph.copy()

    nodes = describe_nodes(graph, neighbors_types=False)
    nodes_to_remove = (
        (nodes['degree'].eq(1) & nodes['type'].eq('DNA'))
        | nodes['subtype'].eq('mRNA')
    )
    nodes_to_remove = nodes.loc[nodes_to_remove, 'node']
    graph.remove_nodes_from(nodes_to_remove)
    return graph


@memory.cache
def build_light_graph(max_workers: int = 2) -> nx.Graph:
    REBUILD_LIGHT_GRAPH = True

    if not REBUILD_LIGHT_GRAPH:
        with (files('biointergraph.static') / "edges_light.tsv.gz").open('rb') as file:
            result = pd.read_csv(file, compression='gzip',
                                 sep='\t', dtype='str')

        result = nx.from_pandas_edgelist(result, edge_attr='dataset')

        return result

    graph = build_main_graph(max_workers=max_workers)
    _lighten_graph(graph, inplace=True)
    graph = _remove_minor_components(graph)
    return graph


def describe_edges(
    graph: nx.Graph, *,
    data: bool = True, explode: bool = False,
    types: bool = False,
    symmetrize: bool = False
) -> pd.DataFrame:
    if data:
        edges = []
        for source, target, attrs in tqdm(graph.edges(data=True), desc='Edges processing: '):
            edges.append((source, target, attrs['dataset']))
        edges = pd.DataFrame(edges, columns=['source', 'target', 'dataset'])
    else:
        edges = pd.DataFrame(graph.edges(), columns=['source', 'target'])

    mask = edges['source'] > edges['target']
    edges.loc[mask, ['source', 'target']] = edges.loc[mask, ['target', 'source']].values
    assert (edges['source'] < edges['target']).all()
    assert edges.shape[0] == graph.number_of_edges()

    if explode:
        assert data
        edges['dataset'] = edges['dataset'].str.split(',')
        edges = edges.explode('dataset')

    if symmetrize:
        edges = pd.concat([
            edges,
            edges.rename(columns={'source': 'target', 'target': 'source'})
        ])
        assert (edges['source'] < edges['target']).mean() == 0.5
        assert edges.shape[0] == 2 * graph.number_of_edges()

    if types:
        edges['source_type'] = _node_id2node_type(edges['source'])
        edges['target_type'] = _node_id2node_type(edges['target'])

    return edges


def node2neighbors(graph: nx.Graph) -> pd.Series:
    result = describe_edges(graph, data=False, symmetrize=True)
    tqdm.pandas()
    result = result.groupby('source')['target'].progress_aggregate(pd.Series.to_list)
    result = result.apply(set)
    return result


def graph2random_walks(graph: nx.Graph, n: int, length: int = 1) -> pd.DataFrame:
    edges = describe_edges(graph, data=False, symmetrize=True)

    result = edges.sample(n, replace=True)
    result = result.rename(columns={'source': 0, 'target': 1})

    tqdm.pandas()
    edges = edges.groupby('source')['target'].progress_aggregate(pd.Series.to_list)

    for i in tqdm(range(1, length)):
        result = result.join(edges.rename(i+1), on=i, how='left', validate='many_to_one')
        tqdm.pandas(leave=False)
        result[i+1] = result[i+1].progress_apply(random.choice)

    return result


def indirect_interactions(graph: nx.Graph, n: int) -> pd.DataFrame:
    result = graph2random_walks(graph, n, 2)[[0,2]]
    result.columns = 'source', 'target'

    result = result[result['source'] != result['target']]
    mask = result['source'] > result['target']
    result.loc[mask, ['source', 'target']] = result.loc[mask, ['target', 'source']].values
    assert (result['source'] < result['target']).all()

    result = result.drop_duplicates()

    neighbors = node2neighbors(graph)
    result['source_neighbors'] = result['source'].map(neighbors)
    result['target_neighbors'] = result['target'].map(neighbors)

    tqdm.pandas()
    result['CN'] = result.progress_apply(
        lambda row: row['source_neighbors'] & row['target_neighbors'],
        axis=1
    )
    result = result.drop(columns=['source_neighbors', 'target_neighbors'])

    result = result.explode('CN')
    result['CN'] = _node_id2node_type(result['CN'])
    result = result.pivot_table(
        index=['source', 'target'], columns='CN',
        aggfunc='size', fill_value=0
    )
    result = result.assign(
        n_types=(result > 0).sum(axis=1),
        entropy=entropy(result, axis=1),
        n_common=result.sum(axis=1)
    )
    result = result.reset_index()

    edges = describe_edges(graph, data=False)
    edges['is_direct'] = True
    result = result.merge(edges, how='left', validate='many_to_one')
    with pd.option_context("future.no_silent_downcasting", True):
        result['is_direct'] = result['is_direct'].fillna(False).infer_objects(copy=False)

    return result


def _merge_singleton_communities(graph: nx.Graph, communities: pd.DataFrame) -> pd.DataFrame:
    member2community = communities.explode('members').set_index('members')
    singleton = member2community[member2community['size'].eq(1)].index

    edges = describe_edges(graph, data=False, symmetrize=True)
    singleton = edges[edges['source'].isin(singleton)]

    singleton = singleton.join(member2community, on='target', how='left', validate='many_to_one')
    singleton = singleton.sort_values('size').drop_duplicates('source', keep='last')

    assert (communities.index == communities['community']).all()
    for _, row in singleton.iterrows():
        communities.loc[row['community'], 'members'].add(row['source'])
        communities.loc[row['community'], 'size'] += 1

    communities = communities[communities['size'].ne(1)].copy()
    assert communities['size'].sum() == graph.number_of_nodes()

    return communities


def _protein_ids2enrichment(ids: list[str]) -> list|float:
    if len(ids) > 10000:
        return float('nan')

    response = requests.post(
        url='https://biit.cs.ut.ee/gprofiler/api/gost/profile/',
        json={
            'organism': 'hsapiens',
            'query': ids,
            'sources': ['GO'],
            'numeric_ns': 'BIOGRID'
        }
    )
    response.raise_for_status()
    response = response.json()
    failed = response['meta']['genes_metadata']['failed']
    if failed:
        print('Failed ids:', *failed)
    result = response['result']
    return result if result else float('nan')


def _community2enrichment(communities: pd.DataFrame) -> pd.Series:
    result = communities.set_index('community', verify_integrity=True)['members'].explode()
    result = result[result.str.startswith('YAPID')]
    result = result.map(yapid2best_id())
    result = result.groupby(level=0).agg(pd.Series.to_list)
    tqdm.pandas()
    result = result.progress_apply(_protein_ids2enrichment)
    return result


def detect_communities(
        graph: nx.Graph,
        merge: bool = False,
        members_types: bool = False,
        enrichment: bool = False,
        nuclear_frac: bool = True,
        **kwargs
    ) -> pd.DataFrame:
    result = nx.community.louvain_communities(graph, **kwargs)
    result = pd.DataFrame({'members': result})
    result = result.reset_index(names='community')

    result['size'] = result['members'].apply(len)
    assert result['size'].sum() == graph.number_of_nodes()

    if merge:
        result = _merge_singleton_communities(graph, result)

    if members_types:
        types = result.explode('members')
        types['type'] = types['members'].str[:5].map({
            'YALID': 'DNA',
            'YAGID': 'RNA',
            'YAPID': 'protein'
        })
        types = types.pivot_table(index='community', columns='type', aggfunc='size', fill_value=0)
        result = pd.concat([result, types], axis=1)
        assert (result['size'] == result[['RNA', 'DNA', 'protein']].sum(axis=1)).all()

    if enrichment:
        result['enrichment'] = _community2enrichment(result)

    if nuclear_frac:
        nuclear_data = result[['community', 'members']].explode('members')
        nuclear_data['nuclear_frac'] = nuclear_data['members'].map(yapid2is_nuclear())
        result = result.join(
            nuclear_data.groupby('community')['nuclear_frac'].mean().fillna(0.0),
            on='community',
            how='left',
            validate='one_to_one'
        )

    return result
