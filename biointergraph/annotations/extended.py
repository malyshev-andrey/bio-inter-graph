import pandas as pd

from ..shared import BED_COLUMNS, memory, GOOGLE_DRIVE_URL, _read_tsv
from ..ids import drop_id_version


def load_extended_annotation(convert2bed: bool = False, **kwargs) -> pd.DataFrame:
    default_kwargs = dict(
        header=None,
        names=[
            'chr' if convert2bed else 'seqid', 'start', 'end',
            'gene_name', 'gene_type', 'strand',
            'source', 'gene_id'
        ],
        usecols=range(8),
        desc='Extended annotation: '
    )
    default_kwargs.update(kwargs)
    result = _read_tsv(
        GOOGLE_DRIVE_URL.format(id='1n2VDbdYe-0di0PVjOKxxk0hZgC914l4e'),
        use_cache=True,
        **default_kwargs
    )

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


@memory.cache
def extended_gene_id2ensembl_gene_id() -> pd.DataFrame:
    result = load_extended_annotation()
    result = result.loc[result['source'].eq('gencode44'), ['gene_id', 'extended_gene_id']]
    result['gene_id'] = drop_id_version(result['gene_id'])
    assert result['gene_id'].str.match(r'^ENSG\d{11}$').all()
    assert not result.duplicated().any()
    return result
