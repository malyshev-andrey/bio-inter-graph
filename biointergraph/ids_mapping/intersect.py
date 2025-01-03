import pandas as pd

from ..shared import memory
from ..ids import drop_id_version
from ..annotations import load_refseq_bed, load_gencode_bed, bed_intersect, load_extended_annotation


@memory.cache
def _load_refseq_data(assembly: str) -> pd.DataFrame:
    refseq_trascripts = load_refseq_bed(assembly=assembly, feature='transcript')
    is_valid = refseq_trascripts['name'].str[:2].isin({'NM', 'NR'})
    print(f'GENCODE/RefSeq intersect: invalid RefSeq IDs frac: {1 - is_valid.mean()}')
    refseq_trascripts = refseq_trascripts[is_valid]

    refseq_data = pd.concat([
        load_refseq_bed(assembly=assembly, feature='gene'),
        refseq_trascripts
    ])
    return refseq_data


def _intersect2pairwise(intersect: pd.DataFrame) -> pd.DataFrame:
    is_proper = intersect['jaccard'] >= 0.8
    print(f'Annotations intersect: improper intersections frac: {1 - is_proper.mean()}')
    intersect = intersect[is_proper]
    intersect['name1'] = drop_id_version(intersect['name1'])
    intersect['name2'] = drop_id_version(intersect['name2'])

    result = intersect[['name1', 'name2']].drop_duplicates()
    return result


@memory.cache
def gencode_refseq_intersect2pairwise(assembly: str) -> pd.DataFrame:
    refseq_data = _load_refseq_data(assembly)
    gencode_data = pd.concat([
        load_gencode_bed(assembly=assembly, feature='gene'),
        load_gencode_bed(assembly=assembly, feature='transcript')
    ])

    result = _intersect2pairwise(bed_intersect(
        refseq_data,
        gencode_data,
        unify_chr_assembly=assembly,
        jaccard=True)
    )
    print(f'GENCODE/RefSeq intersect result for {assembly}: {result.shape}')

    return result


def extended_refseq_intersect2pairwise() -> pd.DataFrame:
    refseq_data = _load_refseq_data('hg38')
    extended_data = load_extended_annotation(convert2bed=True)

    result = _intersect2pairwise(bed_intersect(
        refseq_data,
        extended_data,
        unify_chr_assembly='hg38',
        jaccard=True)
    )
    print(f'Extended/RefSeq intersect result: {result.shape}')
    return result
