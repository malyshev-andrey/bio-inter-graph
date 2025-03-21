import pandas as pd

from ..annotations import load_gencode_bed, load_refseq_bed, bed_intersect
from ..ids_mapping import id2yapid, id2yagid
from ..shared import _read_tsv, memory, BED_COLUMNS


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
    annotation_bed = {
        'gencode': load_gencode_bed,
        'refseq': load_refseq_bed
    }[annotation](assembly='hg38', feature='gene')

    result = bed_intersect(
        result[BED_COLUMNS],
        annotation_bed,
        unify_chr_assembly='hg38',
        strandedness='same',
        jaccard=True,
        how='left'
    )

    no_intersect = result['start2'].eq(-1)
    print(f'POSTAR3 peaks without intersections: {no_intersect.sum()}')
    result = result[~no_intersect]

    peak_id = {f'{c}1': c for c in BED_COLUMNS}
    result = result.rename(columns=peak_id)
    result = result.sort_values('jaccard')
    result = result.drop_duplicates(peak_id.values(), keep='last')

    result['yagid'] = id2yagid(result['name2'])
    assert result['yagid'].str.startswith('YAGID').all()
    result['yapid'] = id2yapid('SYMBOL:' + result['name'])
    assert result['yapid'].str.startswith('YAPID').all()

    result = result[['yagid', 'yapid']]
    result = result.drop_duplicates()

    return result
