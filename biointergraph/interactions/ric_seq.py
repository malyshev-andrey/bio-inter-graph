import numpy as np
import pandas as pd

from ..annotations import load_extended_annotation
from ..shared import memory, CHUNKSIZE, _read_tsv
from ..ids_mapping import id2yagid


def _ricseq_loader(
        id, *,
        chunksize: int|None = CHUNKSIZE,
        pvalue: float|None = None
    ) -> pd.DataFrame:
    url = f'https://drive.usercontent.google.com/download?id={id}&export=download&confirm=t'

    if pvalue is not None:
        filter_func = lambda df: df[df['p_adj'].astype('float') < pvalue]
    else:
        filter_func = lambda df: df

    result = _read_tsv(
        url,
        compression='gzip',
        chunksize=chunksize,
        use_cache=True,
        filter_func=filter_func
    )

    return result


def _load_extended_ricseqlib(**kwargs) -> pd.DataFrame:
    result = _ricseq_loader('1Ie-3g1DQozFyNXvZLwGl4LbKvoBDkm7v', **kwargs)

    result['name'] = result['name'].str.replace('__', ' ')
    gene_id_regex = r'(?:[^ ]+| na)(?: (?:[AB]\.)?\d{1,3})?'
    regex = f'^(?P<gene_id1>{gene_id_regex}) (?P<gene_id2>{gene_id_regex})$'
    result[['gene_id1', 'gene_id2']] = result['name'].str.extract(regex)

    result['gene_id1'] = result['gene_id1'].str.replace(' ', '__')
    result['gene_id2'] = result['gene_id2'].str.replace(' ', '__')

    annotation = load_extended_annotation()
    annotation = annotation.drop('source', axis='columns')

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

    result[['gene_id1', 'gene_id2']] = result[['extended_gene_id1', 'extended_gene_id2']]

    return result


def _load_gencode44_ricseqlib(**kwargs) -> pd.DataFrame:
    result = _ricseq_loader('1zi23ngx_q32zzCV8EaRpCSSGaJKD5x4E', **kwargs)
    result[['gene_id1', 'gene_id2']] = result['name'].str.split('__', regex=False, expand=True)
    return result


def _load_ricpipe(**kwargs) -> pd.DataFrame:
    result = _ricseq_loader('1-2qEi-2EZGpQoLg1povQ0gfSFq33Sh71', **kwargs)
    result[['gene_id1', 'gene_id2']] = result['name'].str.split('_', regex=False, expand=True)
    return result


@memory.cache
def load_ric_seq_data(pvalue: float|None = None, **kwargs) -> pd.DataFrame:
    columns = ['gene_id1', 'gene_id2', 'p_adj']
    if pvalue is not None:
        kwargs['pvalue'] = pvalue

    extended_ricseqlib = _load_extended_ricseqlib(**kwargs)[columns]
    extended_ricseqlib['pipeline'] = 'RICseqlib'
    extended_ricseqlib['annotation'] = 'extended'

    gencode44_ricseqlib = _load_gencode44_ricseqlib(**kwargs)[columns]
    gencode44_ricseqlib['pipeline'] = 'RICseqlib'
    gencode44_ricseqlib['annotation'] = 'gencode44'

    ricpipe = _load_ricpipe(**kwargs)[columns]
    ricpipe['pipeline'] = 'RICpipe'
    ricpipe['annotation'] = 'gencode44'

    result = pd.concat([
        extended_ricseqlib,
        gencode44_ricseqlib,
        ricpipe
    ])

    result['p_adj'] = result['p_adj'].astype('float')

    assert (result['p_adj'] < pvalue).all() or pvalue is None

    print(result.groupby(['pipeline', 'annotation']).size())

    result['yagid1'] = id2yagid(result['gene_id1'])
    result['yagid2'] = id2yagid(result['gene_id2'])
    assert (
        result['yagid1'].str.startswith('YAGID').all() and
        result['yagid2'].str.startswith('YAGID').all()
    )

    result['weight'] = -np.log10(result['p_adj'])

    result = result[['yagid1', 'yagid2', 'weight']]
    result = pd.concat([
        result,
        result.rename(columns={'yagid1': 'yagid2', 'yagid2': 'yagid1'})
    ])
    result = result[result['yagid1'] < result['yagid2']]

    result = result.groupby(['yagid1', 'yagid2'], as_index=False, observed=True)['weight'].max()

    return result
