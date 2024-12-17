import pandas as pd
from tqdm.auto import tqdm

from ..annotations import load_extended_annotation


def load_extended_ric_seq_data(
        chunksize: int|None = None,
        pvalue: float = 0.2
    ) -> pd.DataFrame:
    url = 'https://drive.usercontent.google.com/download?id=1mx0LHzmZ8Zfg-wQWJR1Ok1JztwOUZqZW&export=download&confirm=t'

    default_kwargs = dict(
        sep='\t',
        compression='zip',
        dtype='str'
    )
    if chunksize is None:
        result = pd.read_csv(url, **default_kwargs)
    else:
        result = []
        with tqdm(desc=url) as progress_bar:
            for chunk in pd.read_csv(url, chunksize=chunksize, **default_kwargs):
                progress_bar.update(chunk.shape[0])
                result.append(chunk)
        result = pd.concat(result)

    result = result[result['p_adj'].astype('float') < pvalue].copy()

    result['name'] = result['name'].str.replace('__', ' ')
    gene_id_regex = r'(?:[^ ]+| na)(?: (?:[AB]\.)?\d{1,3})?'
    regex = f'^(?P<gene_id1>{gene_id_regex}) (?P<gene_id2>{gene_id_regex})$'
    result[['gene_id1', 'gene_id2']] = result['name'].str.extract(regex)

    result['gene_id1'] = result['gene_id1'].str.replace(' ', '__')
    result['gene_id2'] = result['gene_id2'].str.replace(' ', '__')

    annotation = load_extended_annotation(
        names=[
            'seqid', 'start', 'end', 'gene_name',
            'gene_type', 'strand', 'source', 'gene_id'
        ]
    ).drop('source', axis='columns')

    map1 = {c: f'{c}1' for c in annotation.columns}
    map2 = {c: f'{c}2' for c in annotation.columns}

    result = result.merge(
        annotation.rename(columns=map1),
        on=[c for c in map1.values() if c != 'extended_gene_id1'],
        how='left',
        validate='many_to_one'
    )
    assert not result['extended_gene_id1'].isna().any()

    result = result.merge(
        annotation.rename(columns=map2),
        on=[c for c in map2.values() if c != 'extended_gene_id2'],
        how='left',
        validate='many_to_one'
    )
    assert not result['extended_gene_id2'].isna().any()

    return result
