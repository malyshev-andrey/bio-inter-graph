import pandas as pd

from ..annotations import unify_chr
from ..shared import UNIFY_BIOTYPES, memory


@memory.cache
def entrezgene_id_info(**kwargs) -> pd.DataFrame:
    url = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz'
    result = pd.read_csv(url, sep='\t', dtype='str', **kwargs)
    result = result[result['#tax_id'].eq('9606')]
    assert result['GeneID'].str.isdigit().all()

    columns = {
        'GeneID': 'entrezgene_id',
        'chromosome': 'chr',
        'type_of_gene': 'biotype'
    }
    result = result.rename(columns=columns)
    result = result[columns.values()]

    result['chr'] = unify_chr(result['chr'], assembly='hg38')
    assert result['chr'].eq('chrM').any() and not result['chr'].eq('chrMT').any()

    assert result['entrezgene_id'].is_unique

    return result


def entrezgene_id2biotype(ids: pd.Series|None = None) -> pd.Series:
    result = entrezgene_id_info()

    result = result.set_index('entrezgene_id', verify_integrity=True)
    result = result['biotype']

    invalid_biotypes = ['biological-region', 'ncRNA', 'unknown', 'other']
    result = result.replace(invalid_biotypes, float('nan'))
    result = result.replace(UNIFY_BIOTYPES)

    if ids is not None:
        result = ids.map(result)
    return result
