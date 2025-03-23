from concurrent.futures import ThreadPoolExecutor, as_completed
from time import sleep
import importlib.resources

import requests
import pandas as pd
import networkx as nx
from tqdm.auto import tqdm

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
from ..ids_mapping.protein import _build_yapid_graph


def _remove_minor_components(graph):
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
            result = pd.read_csv(
                file, compression='gzip',
                sep='\t', header=None,
                names=['source', 'target']
            )

        result = nx.from_pandas_edgelist(result)

        return result


    data = [
        (load_encode_eclip_data, dict(assembly='hg38', annotation='gencode', cell_line='K562')),
        (load_encode_rip_data, dict(annotation='gencode', cell_line='K562')),
        (load_encode_iclip_data, dict(annotation='gencode', cell_line='K562')),
        (load_postar3_data, dict(species='human', cell_line='K562', annotation='gencode')),
        (load_frip_seq_data, dict()),
        (load_ric_seq_data, dict(pvalue=0.05)),
        (load_karr_seq_data, dict(pvalue=0.05)),
        (load_intact_interactions, dict()),
        (load_biogrid_interactions, dict()),
        (load_string_interactions, dict(min_score=700)),
        (load_encode_chip_seq_data, dict(assembly='hg38', cell_line='K562')),
        (load_redc_redchip_data, dict()),
        (load_gtrd_chip_seq_data, dict(cell_line='K562'))
    ]

    tqdm_kwargs = dict(total=len(data), unit='source', desc='Collecting data: ')
    if max_workers > 1:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = []
            for func, kwargs in data:
                futures.append(executor.submit(func, **kwargs))
                sleep(1)

            data = [f.result() for f in tqdm(as_completed(futures), **tqdm_kwargs)]

    else:
        data = [func(**kwargs) for func, kwargs in tqdm(data, **tqdm_kwargs)]

    for df in data:
        assert df.shape[1] == 2
        df.columns = 'source', 'target'

    data = pd.concat(data)

    assert data.shape[1] == 2
    regex = r'^YA[LPG]ID\d{7}$'
    assert data['source'].str.match(regex).all()
    assert data['target'].str.match(regex).all()

    graph = nx.from_pandas_edgelist(data)

    graph = _remove_minor_components(graph)

    return graph


def describe_graph(graph):
    print(f'Nodes: {graph.number_of_nodes()}')
    print(f'Edges: {graph.number_of_edges()}')
    print('Degrees distribution:')
    stats = pd.Series(dict(graph.degree())).describe()
    print('\t' + str(stats[:10]).replace('\n', '\n\t'))
    diameter_sample = [nx.approximation.diameter(graph) for _ in range(10)]
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

    return result


def describe_nodes(graph) -> pd.DataFrame:
    result = pd.DataFrame(graph.degree(), columns=['node', 'degree'])
    assert result['node'].str.match(r'^YA[GP]ID\d{7}$').all()
    result['is_protein'] = result['node'].str.startswith('YAPID')
    return result
