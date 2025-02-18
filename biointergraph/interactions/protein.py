import pandas as pd
from tqdm.auto import tqdm


from ..shared import CHUNKSIZE, memory


@memory.cache
def ensembl_protein2biogrid_id(
        chunksize: int | None = CHUNKSIZE
    ) -> pd.DataFrame:
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
    assert not result.duplicated().any()
    result = result.set_index('IDENTIFIER_VALUE', verify_integrity=True)
    result = result['BIOGRID_ID']
    return result
