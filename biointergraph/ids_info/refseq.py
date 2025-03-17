import pandas as pd

from ..shared import UNIFY_BIOTYPES, memory
from ..annotations import read_feature_table, unify_chr
from .entrez import entrezgene_id2biotype


def expand_attributes(ft: pd.DataFrame) -> pd.DataFrame:
    n = ft.shape[0]
    result = ft['attributes'].str.split(';')
    result = result.explode()
    result = result.str.split('=', expand=True)
    result = result.pivot(columns=0, values=1)
    result = pd.concat([ft, result], axis='columns')
    assert result.shape[0] == n
    return result


@memory.cache
def refseq_transcript_id_info(**kwargs) -> pd.DataFrame:
    global CACHE
    url = 'https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/historical/GRCh38/GCF_000001405.40-RS_2024_08_historical/GCF_000001405.40-RS_2024_08_genomic.gff.gz'
    type_filter = lambda df: df[~df['type'].isin({
        'gene', 'pseudogene',
        'exon', 'CDS'
    })]
    result = read_feature_table(url, filter_func=type_filter) if CACHE is None else CACHE
    CACHE = result

    result = expand_attributes(result)
    regex = r'^rna-(?P<accession>N[MR]_\d+)\.(?P<version>\d+)(?:-(?P<subversion>\d+))?$'
    result = result[result['ID'].str.match(regex)]
    result = pd.concat([result, result['ID'].str.extract(regex)], axis='columns')
    result['version'] = result['version'].astype('int')
    assert not result['subversion'].eq('0').any()
    result['subversion'] = result['subversion'].fillna('0').astype('int')

    result = result.sort_values(['version', 'subversion'], ascending=False)
    result = result.drop_duplicates('accession', keep='first')

    gene_id_regex = r'GeneID:(\d+)'
    assert result['Dbxref'].str.count(gene_id_regex).eq(1).all()
    result['gene_id'] = result['Dbxref'].str.extract(gene_id_regex, expand=False)
    result['gene_type'] = entrezgene_id2biotype(result['gene_id'])

    columns = {
        'chr': 'chr', 'start': 'start', 'end': 'end', 'strand': 'strand',
        'type': 'type', 'gbkey': 'gbkey', 'pseudo': 'pseudo', 'gene_type': 'gene_type',
        'product': 'product',
        'accession': 'refseq_transcript_id'
    }
    result = result.rename(columns=columns)
    result = result[columns.values()]

    result['chr'] = unify_chr(result['chr'], assembly='hg38')

    assert result['refseq_transcript_id'].is_unique

    return result


def refseq_transcript_id2biotype(ids: pd.Series|None = None) -> pd.Series:
    result = refseq_transcript_id_info()

    result['biotype'] = result['gene_type'].where(
        result['gbkey'].eq('misc_RNA') & result['type'].eq('transcript'),
        result['type']
    )
    result['biotype'] = result['biotype'].where(
        ~result['type'].isin({'primary_transcript', 'miRNA'}),
        'miRNA'
    )
    assert result.loc[result['biotype'].eq('miRNA'), 'product'].str.startswith('microRNA ').all()
    result['biotype'] = result['biotype'].where(
        ~(result['pseudo'].eq('true') | result['gene_type'].eq('pseudogene')),
        'pseudogene'
    )

    result = result.set_index('refseq_transcript_id', verify_integrity=True)
    result = result['biotype']
    result = result.replace(UNIFY_BIOTYPES)

    if ids is not None:
        result = ids.map(result)
    return result
