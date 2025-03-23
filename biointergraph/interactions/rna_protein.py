import pandas as pd

from ..annotations import load_gencode_bed, load_refseq_bed, bed_intersect
from ..ids_mapping import id2yapid, id2yagid
from ..shared import _read_tsv, memory, BED_COLUMNS
from .main import _annotate_peaks


def _load_postar3_peaks(species: str = 'human', **kwargs) -> pd.DataFrame:
    default_kwargs = dict(
        header=None,
        names=[
            'chr', 'start', 'end', 'peak_id', 'strand', 'name',
            'method', 'cell_line', 'accessions', 'score'
        ]
    )
    default_kwargs.update(kwargs)

    result = _read_tsv(
        f'https://cloud.tsinghua.edu.cn/d/8133e49661e24ef7a915/files/?dl=1&p={species}.txt.gz',
        **default_kwargs
    )
    return result


@memory.cache
def load_postar3_data(species: str, cell_line: str, annotation: str, **kwargs) -> pd.DataFrame:
    result = _load_postar3_peaks(
        species=species,
        filter_func=lambda df: df[df['cell_line'].eq(cell_line)],
        **kwargs
    )
    annotation = {
        'gencode': load_gencode_bed,
        'refseq': load_refseq_bed
    }[annotation](assembly='hg38', feature='gene')

    result = _annotate_peaks(
        result[BED_COLUMNS].copy(), annotation,
        assembly='hg38',
        desc='POSTAR3',
        convert_ids=True
    )

    return result

c
def load_frip_seq_data() -> pd.DataFrame:
    result = pd.concat([
        pd.read_csv(
            'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE67nnn/GSE67963/suppl/GSE67963_significance_calls_genes.txt.gz',
            sep='\t'
        ),
        pd.read_csv(
            'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE67nnn/GSE67963/suppl/GSE67963_significance_calls_isoforms.txt.gz',
            sep='\t'
        )
    ])
    result = result[result['significant'].eq('yes')]
    result = pd.concat([
        result[
            result['sample_1'].eq('input') &
            (result['value_1'] < result['value_2'])
        ],
        result[
            result['sample_2'].eq('input') &
            (result['value_1'] > result['value_2'])
        ]
    ])

    result['sample_1'] = result['sample_2'].where(
        result['sample_1'].eq('input'),
        result['sample_1']
    )
    assert not result['sample_1'].eq('input').any()

    result['gene_id'] = result['gene_id'].str.split(',')
    result = result.explode('gene_id')
    result['yagid'] = id2yagid(result['gene_id'])
    result = result[result['yagid'].str.startswith('YAGID')]

    result['sample_1'] = result['sample_1'].str.upper()
    result['sample_1'] = result['sample_1'].replace({
        'CBP': 'CREBBP',
        'CBX': 'CBX3',
        'HNRNPH': 'HNRNPH1',
        'IMP1': 'IGF2BP1',
        'HUR': 'ELAVL1',
        'LSD1': 'KDM1A',
        'PCAF': 'KAT2B',
        'PABP': 'PABPC1'
    })
    result['yapid'] = id2yapid('SYMBOL:' + result['sample_1'])
    result = result[result['yapid'].str.startswith('YAPID')]

    result = result[['yagid', 'yapid']]
    result = result.drop_duplicates()

    return result
