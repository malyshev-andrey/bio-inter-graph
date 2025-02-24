import pandas as pd


from ..shared import memory
from ..ids_mapping import id2yapid


@memory.cache
def load_biogrid_interactions() -> pd.DataFrame:
    url = 'https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-MV-Physical-LATEST.tab3.zip'
    result = pd.read_csv(url, sep='\t', dtype='str')
    result = result[
        result['Organism ID Interactor A'].eq('9606') &
        result['Organism ID Interactor B'].eq('9606') &
        ~result['Experimental System'].isin(['Protein-RNA', 'Affinity Capture-RNA'])
    ]
    assert (
        result['Organism Name Interactor A'].eq('Homo sapiens') &
        result['Organism Name Interactor B'].eq('Homo sapiens')
    ).all()

    result['yapid1'] = id2yapid(result['BioGRID ID Interactor A'])
    result['yapid2'] = id2yapid(result['BioGRID ID Interactor B'])

    ids = ['yapid1', 'yapid2']
    result = result[ids]
    swap = dict(zip(ids, ids[::-1]))
    result = pd.concat([
        result,
        result.rename(columns=swap)
    ])
    result = result[result['yapid1'] < result['yapid2']]
    result = result.drop_duplicates()
    return result


@memory.cache
def load_string_interactions() -> pd.DataFrame:
    url = 'https://stringdb-downloads.org/download/stream/protein.physical.links.v12.0/9606.protein.physical.links.v12.0.onlyAB.tsv.gz'
    result = pd.read_csv(url, sep='\t')
    assert (
        result['protein1'].str.match('^9606.ENSP\d{11}$') &
        result['protein2'].str.match('^9606.ENSP\d{11}$')
    ).all()
    result['protein1'] = result['protein1'].str.removeprefix('9606.')
    result['protein2'] = result['protein2'].str.removeprefix('9606.')

    result = result[['protein1', 'protein2']]
    result = result.rename(columns={
        'protein1': 'ensembl_protein1',
        'protein2': 'ensembl_protein2'
    })
    assert not result.duplicated().any()
    return result
