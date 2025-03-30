import pandas as pd
import networkx as nx
from tqdm.auto import tqdm

from .graph import describe_nodes
from .main import summarize_pairwise


def graph2rna_protein(graph: nx.Graph) -> pd.DataFrame:
    result = describe_nodes(graph, subtypes=False, neighbors_types=False)
    result = result.loc[result['type'].eq('DNA'), ['node', 'neighbors']].copy()
    n_nodes = result.shape[0]
    tqdm.pandas(desc='Neighbors processing: ')
    result['RNA'] = result['neighbors'].progress_apply(
        lambda nodes: [n for n in nodes if n.startswith('YAGID')]
    )
    result['protein'] = result['neighbors'].progress_apply(
        lambda nodes: [n for n in nodes if n.startswith('YAPID')]
    )
    result = result.drop(columns='neighbors')
    result = result.explode('RNA').explode('protein')
    assert result['node'].nunique() == n_nodes

    result = result[~result.isna().any(axis=1)]
    assert not result.duplicated().any()

    result = summarize_pairwise(result, ['RNA', 'protein'], nunique=('node', 'nunique'))
    assert (result['nunique'] == result['size']).all()

    edges = pd.DataFrame(graph.edges(), columns=['source', 'target'])
    swap_mask = edges['source'] > edges['target']
    edges.loc[swap_mask, ['source', 'target']] = edges.loc[swap_mask, ['target', 'source']].values
    assert (edges['source'] < edges['target']).all()

    assert not (
        edges['source'].str.startswith('YAPID') &
        edges['target'].str.startswith('YAGID')
    ).any()

    edges = edges[
        edges['source'].str.startswith('YAGID') &
        edges['target'].str.startswith('YAPID')
    ]
    edges = edges.rename(columns={'source': 'RNA', 'target': 'protein'})
    edges['is_direct'] = True

    result = result.merge(edges, how='left', validate='one_to_one')
    result['is_direct'] = result['is_direct'].fillna(False).infer_objects()

    return result
