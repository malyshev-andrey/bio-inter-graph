import importlib.resources

import pandas as pd
import networkx as nx

from ..shared import memory, _read_tsv


def _retrieve_string_ids() -> pd.Series:
    url = 'https://stringdb-downloads.org/download/stream/protein.physical.links.v12.0/9606.protein.physical.links.v12.0.onlyAB.tsv.gz'
    result = pd.read_csv(url, sep='\t')
    result = pd.concat([
        result['protein1'],
        result['protein2']
    ])
    assert result.str.match(r'^9606\.ENSP\d{11}$').all()
    result = result.str.removeprefix('9606.')
    result = result.drop_duplicates()
    return result


def _retrieve_intact_ids() -> pd.Series:
    result = _read_tsv(
        'https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/species/human.txt',
        usecols=['#ID(s) interactor A', 'ID(s) interactor B']
    )
    result = pd.concat([
        result['#ID(s) interactor A'],
        result['ID(s) interactor B']
    ])
    result = result[result.str.startswith('uniprotkb:')]
    result = result.str.removeprefix('uniprotkb:')
    result = result.str.split('-', expand=True)[0]
    regex = r'^([A-Z0-9]{6}|[A-Z0-9]{10})$'
    assert result.str.match(regex).all()
    return result


def _string_mapping() -> pd.DataFrame:
    result = _read_tsv(
        'https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz',
        filter_func=lambda df: df[
            df['source'].isin({'UniProt_AC', 'UniProt_GN_Name'})
        ]
    )

    result['alias'] = result['alias'].where(
        ~result['source'].eq('UniProt_GN_Name'),
        'SYMBOL:' + result['alias']
    )

    regex = r'^(SYMBOL:.+|[A-Z0-9]{6}|[A-Z0-9]{10})$'
    assert result['alias'].str.match(regex).all()

    regex = r'^9606\.ENSP\d{11}$'
    assert result['#string_protein_id'].str.match(regex).all()
    result['#string_protein_id'] = result['#string_protein_id'].str.removeprefix('9606.')

    result = result[['#string_protein_id', 'alias']]
    result.columns = 'source', 'target'

    return result


def _biogrid_mapping() -> pd.DataFrame:
    result = _read_tsv(
        'https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-IDENTIFIERS-LATEST.tab.zip',
        filter_func=lambda df: df[
            df['ORGANISM_OFFICIAL_NAME'].eq('Homo sapiens') &
            df['IDENTIFIER_TYPE'].isin({
                'TREMBL',
                'UNIPROT-ACCESSION',
                'SWISS-PROT',
                'ENSEMBL PROTEIN',
                'OFFICIAL SYMBOL'
            })
        ],
        skiprows=27
    )

    assert result['BIOGRID_ID'].str.isdigit().all()
    assert not result['IDENTIFIER_VALUE'].str.isdigit().any()

    result['IDENTIFIER_VALUE'] = result['IDENTIFIER_VALUE'].where(
        ~result['IDENTIFIER_TYPE'].eq('OFFICIAL SYMBOL'),
        'SYMBOL:' + result['IDENTIFIER_VALUE']
    )
    regex = r'^(ENSP\d{11}|SYMBOL:.+|[A-Z0-9]{6}|[A-Z0-9]{10})$'
    assert result['IDENTIFIER_VALUE'].str.match(regex).all()

    result = result[['BIOGRID_ID', 'IDENTIFIER_VALUE']]
    result.columns = 'source', 'target'

    return result


@memory.cache
def _build_yapid_graph():
    REBUILD_YAPID_MAPPING = False

    if not REBUILD_YAPID_MAPPING:
        with importlib.resources.open_text('bio-inter-graph.static', 'id2yapid.json') as file:
            result = pd.read_json(file, typ='series')
        return result

    mapping = pd.concat([_string_mapping(), _biogrid_mapping()])
    assert mapping.shape[1] == 2

    yapid_graph = nx.from_pandas_edgelist(mapping)
    yapid_graph.add_nodes_from(_retrieve_string_ids())
    yapid_graph.add_nodes_from(_retrieve_intact_ids())

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


def id2yapid(ids: pd.Series|None = None) -> pd.Series:
    result = _build_yapid_graph()
    if ids is not None:
        result = ids.map(result).combine_first(ids)
    return result
