import pkg_resources

import numpy as np
import pandas as pd
import networkx as nx

from ..shared import REBUILD_YAPID_MAPPING

def _build_yapid_graph():
    if not REBUILD_YAPID_MAPPING:
        filepath = pkg_resources.resource_filename(
            'bio-inter-graph',
            'static/id2yapid.json'
        )
        result = pd.read_json(filepath, typ='series')
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

    connected_components = nx.connected_components(yapid_graph)

    result = {}
    for component in connected_components:
        biogrid_id = float('-inf')
        for node in component:
            if node.isdigit():
                biogrid_id = max(biogrid_id, int(node))
        assert float('-inf') < biogrid_id < 1e7

        for node in component:
            result[node] = 'YAPID' + str(biogrid_id).zfill(7)
    result = pd.Series(result)

    return result


def id2yapid(ids: pd.Series) -> pd.Series:
    mapping = _build_yapid_graph()
    ids = ids.map(mapping).combine_first(ids)
    return ids
