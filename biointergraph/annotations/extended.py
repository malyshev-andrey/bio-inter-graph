import pandas as pd

from ..shared import BED_COLUMNS


def load_extended_annotation(convert2bed: bool = False, **kwargs) -> pd.DataFrame:
    url = 'https://drive.usercontent.google.com/download?id=1n2VDbdYe-0di0PVjOKxxk0hZgC914l4e&export=download&confirm=t'

    default_kwargs = dict(
        sep='\t',
        header=None,
        names=[
            'chr' if convert2bed else 'seqid', 'start', 'end',
            'gene_name', 'gene_type', 'strand',
            'source', 'gene_id'
        ],
        usecols=range(8),
        dtype='str'
    )
    default_kwargs.update(kwargs)

    result = pd.read_csv(url, **default_kwargs)

    result['extended_gene_id'] = result.index.astype('str')
    result['extended_gene_id'] = result['extended_gene_id'].str.zfill(7)
    result['extended_gene_id'] = 'EXTG' + result['extended_gene_id']

    assert result['start'].str.isdigit().all()
    assert result['end'].str.isdigit().all()
    assert result['strand'].isin({'+', '-'}).all()
    assert result['extended_gene_id'].is_unique

    assert (result['start'].astype('int') < result['end'].astype('int')).all()

    if not convert2bed:
        return result

    result['start'] = result['start'].astype('int')
    result['end'] = result['end'].astype('int')
    result['name'] = result['extended_gene_id']
    result['score'] = '.'

    return result[BED_COLUMNS]
