import pandas as pd
from tqdm.auto import tqdm


from ..shared import CHUNKSIZE, memory


@memory.cache
def ensembl_protein2biogrid_id(
        chunksize: int | None = CHUNKSIZE
    ) -> pd.Series:
    url = 'https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-IDENTIFIERS-LATEST.tab.zip'
    filter_func = lambda df: df[
        df['ORGANISM_OFFICIAL_NAME'].eq('Homo sapiens') &
        df['IDENTIFIER_TYPE'].eq('ENSEMBL PROTEIN')
    ]

    read_csv_kwargs = dict(sep='\t', skiprows=27, dtype='str')

    if chunksize is None:
        result = filter_func(pd.read_csv(url, **read_csv_kwargs))
    else:
        result = []
        with tqdm() as progress_bar:
            for chunk in pd.read_csv(url, chunksize=chunksize, **read_csv_kwargs):
                progress_bar.update(chunk.shape[0])
                result.append(filter_func(chunk))
        result = pd.concat(result)

    result = result[['BIOGRID_ID', 'IDENTIFIER_VALUE']]
    result = result.rename(columns={
        'BIOGRID_ID': 'biogrid_id',
        'IDENTIFIER_VALUE': 'ensembl_protein'
    })
    assert not result.duplicated().any()
    result = result.set_index('ensembl_protein', verify_integrity=True)
    result = result['biogrid_id']
    return result


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

    result = result[['BioGRID ID Interactor A', 'BioGRID ID Interactor B']]
    result = result.rename(columns={
        'BioGRID ID Interactor A': 'biogrid_id1',
        'BioGRID ID Interactor B': 'biogrid_id2'
    })
    result = result.drop_duplicates()
    return result


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
