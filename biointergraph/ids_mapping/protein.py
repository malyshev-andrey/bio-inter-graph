import importlib.resources

import numpy as np
import pandas as pd
import networkx as nx

from ..shared import memory


def _retrieve_string_ids() -> pd.Series:
    url = 'https://stringdb-downloads.org/download/stream/protein.physical.links.v12.0/9606.protein.physical.links.v12.0.onlyAB.tsv.gz'
    result = pd.read_csv(url, sep='\t')
    result = pd.concat([
        result['protein1'],
        result['protein2']
    ])
    result = result.str.removeprefix('9606.')
    result = result.drop_duplicates()
    return result


@memory.cache
def _build_yapid_graph():
    REBUILD_YAPID_MAPPING = True

    if not REBUILD_YAPID_MAPPING:
        with importlib.resources.open_text('bio-inter-graph.static', 'id2yapid.json') as file:
            result = pd.read_json(file, typ='series')
        return result

    read_csv_kwargs = dict(sep='\t', skiprows=27, dtype='str')
    url = 'https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-IDENTIFIERS-LATEST.tab.zip'
    filter_func = lambda df: df[
        df['ORGANISM_OFFICIAL_NAME'].eq('Homo sapiens') &
        df['IDENTIFIER_TYPE'].isin([
            'TREMBL',
            'UNIPROT-ACCESSION',
            'SWISS-PROT',
            'ENSEMBL PROTEIN',
            'OFFICIAL SYMBOL'
        ])
    ]

    mapping = filter_func(pd.read_csv(url, **read_csv_kwargs))

    assert mapping['BIOGRID_ID'].str.isdigit().all()
    assert not mapping['IDENTIFIER_VALUE'].str.isdigit().any()

    mapping['IDENTIFIER_VALUE'] = np.where(
        mapping['IDENTIFIER_TYPE'].eq('OFFICIAL SYMBOL'),
        'SYMBOL:' + mapping['IDENTIFIER_VALUE'],
        mapping['IDENTIFIER_VALUE']
    )

    mapping = mapping[['BIOGRID_ID', 'IDENTIFIER_VALUE']]
    mapping.columns = ['source', 'target']

    yapid_graph = nx.from_pandas_edgelist(mapping)
    yapid_graph.add_nodes_from(_retrieve_string_ids())

    connected_components = nx.connected_components(yapid_graph)

    result = {}
    for i, component in enumerate(connected_components):
        assert i < 1e7
        for node in component:
            result[node] = 'YAPID' + str(i).zfill(7)
    result = pd.Series(result)

    n_components = result.nunique()
    print(f'Build YAPID graph: {len(result)} ids, {n_components} components')

    return result


def id2yapid(ids: pd.Series) -> pd.Series:
    mapping = _build_yapid_graph()
    ids = ids.map(mapping).combine_first(ids)
    return ids
