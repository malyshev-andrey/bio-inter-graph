import pandas as pd

from ..shared import GOOGLE_DRIVE_URL, memory
from ..annotations import load_extended_annotation, bed_intersect, load_ChromHMM_annotation
from ..ids_mapping import id2yagid


@memory.cache
def load_rna_chrom_data(**kwargs) -> pd.DataFrame:
    url = GOOGLE_DRIVE_URL.format(id='1nkg0Iofz8azz6BEfISWXG_DNQMlHEbi6')
    default_kwargs = dict(sep='\t', compression='zip', dtype='str')
    default_kwargs.update(kwargs)
    result = pd.read_csv(url, **default_kwargs)

    names_map = {
        'seqid': 'chr1', 'start': 'gene_start',
        'end': 'gene_end', 'strand': 'strand1',
        'extended_gene_id': 'name'
    }
    annotation = load_extended_annotation()
    annotation = annotation.rename(columns=names_map)

    result = result.merge(annotation, how='left', validate='many_to_one')
    assert not result['name'].isna().any()

    cols = ['chr', 'start', 'end', 'name']
    result = result.rename(columns={f'{c}2': c for c in cols})
    result = result[cols]

    result = bed_intersect(
        result, load_ChromHMM_annotation()[cols],
        strandedness = None,
        unify_chr_assembly = 'hg38',
        jaccard = True,
        how='left'
    )
    no_intersections = result['name2'].eq('-1')
    print(f'Peaks without intersections: {no_intersections.sum()}')
    result = result[~no_intersections]

    peak_id = ['chr', 'start1', 'end1', 'name1']
    result = result.sort_values('jaccard', ascending=False)
    result = result.drop_duplicates(peak_id, keep='first')

    result['yagid'] = id2yagid(result['name1'])
    assert result['yagid'].str.startswith('YAGID').all()

    result = result.rename(columns={'name2': 'yalid'})
    result = result[['yagid', 'yalid']].drop_duplicates()

    return result
