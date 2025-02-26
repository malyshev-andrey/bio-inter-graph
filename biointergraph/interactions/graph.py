import pandas as pd
import networkx as nx

from .ENCODE import encode_eCLIP2pairwise
from .karr_seq import load_karr_seq_data
from .ric_seq import load_ric_seq_data
from .protein import load_IntAct_interactions, load_biogrid_interactions, load_string_interactions
from ..shared import memory


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
def build_main_graph():
    data = [
        encode_eCLIP2pairwise(assembly='hg38', annotation='gencode', pvalue=0.05),
        load_ric_seq_data(0.05),
        load_karr_seq_data(0.05),
        load_IntAct_interactions(),
        load_biogrid_interactions(),
        load_string_interactions()
    ]

    for df in data:
        assert df.shape[1] == 2
        df.columns = 'source', 'target'

    data = pd.concat(data)
    assert data.shape[1] == 2
    regex = r'^YA[PG]ID\d{7}$'
    assert (
        data['source'].str.match(regex).all() and
        data['target'].str.match(regex).all()
    )

    graph = nx.from_pandas_edgelist(data)

    graph = _remove_minor_components(graph)

    return graph


def describe_graph(graph):
    print(f'Nodes: {graph.number_of_nodes()}')
    print(f'Edges: {graph.number_of_edges()}')
    print('Degrees distribution:')
    print(pd.Series(dict(graph.degree())).describe())
    diameter_sample = [nx.approximation.diameter(graph) for _ in range(10)]
    print(f'Diameter: {min(diameter_sample)}-{max(diameter_sample)}')


def detect_communities(graph) -> pd.DataFrame:
    communities = nx.community.louvain_communities(graph, resolution=2)
    result = []
    for community in communities:
        community = list(community)
        result.append({
            'nodes': community,
            'size': len(community),
            'protein_frac': pd.Series(community).str.startswith('YAPID').mean()
        })
    result = pd.DataFrame(result)

    return result


def describe_nodes(graph) -> pd.DataFrame:
    result = pd.DataFrame(graph.degree(), columns=['node', 'degree'])
    assert result['node'].str.match(r'^YA[GP]ID\d{7}$').all()
    result['is_protein'] = result['node'].str.startswith('YAPID')
    return result
