import pandas as pd
import networkx as nx
from tqdm.auto import tqdm

from .graph import describe_nodes, describe_edges
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


def graph_datasets_stats(graph: nx.Graph) -> pd.DataFrame:
    edges = describe_edges(graph, explode=True, types=True)
    assert (edges['source'] < edges['target']).all()

    result = edges.groupby(
        ['source_type', 'target_type', 'dataset'],
        as_index=False,
        observed=True
    ).agg(
        interactions=('dataset', 'size'),
        n_sources=('source', 'nunique'),
        n_targets=('target', 'nunique')
    )
    index =                ['Source',     'Protocol/Database',       'Cell line', 'Annotation',        'Assembly']
    metadata = {
        'ENCODE ChIP-seq': ['ENCODE',     'ChIP-seq',                'K562',      'GENCODE',           'hg38'    ],
        'Red-C & RedChIP': ['GEO/BaRDIC', 'Red-C & RedChIP (input)', 'K562',      'Extended/ChromHMM', 'hg38'    ],
        'GTRD':            ['GTRD',       'ChIP-seq meta-clusters',  'K562',      'ChromHMM',          'hg38'    ],
        'RIC-seq':         ['GSE190214',  'RIC-seq',                 'K562',      'Extended/GENCODE',  'hg38'    ],
        'POSTAR3':         ['POSTAR3',    'Database',                'K562',      'GENCODE',           'hg38'    ],
        'KARR-seq':        ['GSE166155',  'KARR-seq',                'K562',      'RefSeq',            'hg19'    ],
        'ENCODE eCLIP':    ['ENCODE',     'eCLIP',                   'K562',      'GENCODE',           'hg38'    ],
        'IntAct':          ['IntAct',     'Database',                'Human',     '-',                 '-'       ],
        'BioGRID':         ['BioGRID',    'Database',                'Human',     '-',                 '-'       ],
        'STRING':          ['STRING',     'Database',                'Human',     '-',                 '-'       ],
        'fRIP-seq':        ['GSE67963',   'fRIP-seq',                'K562',      'GENCODE',           'hg19'    ],
        'ENCODE RIP':      ['ENCODE',     'RIP-seq & RIP-chip',      'K562',      'GENCODE',           'hg19'    ],
        'ENCODE iCLIP':    ['ENCODE',     'iCLIP',                   'K562',      'GENCODE',           'hg19'    ]
    }
    metadata = pd.DataFrame(metadata, index=index).T

    result = result.join(metadata, on='dataset', how='left', validate='one_to_one')
    result = result.drop(columns='dataset')
    result = result.sort_values('interactions', ascending=False)
    return result


def graph_datasets_matrix(graph: nx.Graph) -> pd.DataFrame:
    edges = describe_edges(graph, explode=True)

    edges['_temp'] = True
    result = edges.pivot(
        index=['source', 'target'],
        columns='dataset',
        values='_temp'
    ).reset_index(drop=True).fillna(False)

    return result
