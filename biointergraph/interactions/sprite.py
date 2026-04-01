import gzip
import re
from pathlib import Path

import pandas as pd
from tqdm.auto import tqdm

from ..annotations import (
    best_left_intersect,
    load_chromhmm_annotation,
    load_extended_annotation
)
from ..ids_mapping import id2yagid
from ..shared import BED_COLUMNS, memory, remote_file2local


RDSPRITE_DATA_URL = (
    'https://www.ncbi.nlm.nih.gov/geo/download/'
    '?acc=GSM8173720&format=file&file=GSM8173720_H1.RDSPRITE.comboall.txt.gz'
)

DNA_REGEX = re.compile(r'^DPM\[[+-]\]_(chr[^:]+):(\d+)-(\d+)$')
RNA_COORD_REGEX = re.compile(r'_(chr[^:]+):(\d+)-(\d+)$')
RNA_REFSEQ_REGEX = re.compile(r'(N[MR]_\d+(?:\.\d+)?)')


def _make_unique(values: list[str]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(values))


def _as_local_path(url: str) -> Path:
    local_url = remote_file2local(url)
    assert local_url.startswith('file://')
    return Path(local_url.removeprefix('file://'))


def _parse_rdsprite_clusters(filepath: Path) -> pd.DataFrame:
    result = []

    with gzip.open(filepath, 'rt') as file:
        for cluster_index, line in enumerate(
                tqdm(file, desc='READING RD-SPRITE clusters', unit='cluster')
            ):
            fields = line.rstrip('\n').split('\t')
            cluster_id, members = fields[0], fields[1:]

            rna_members = [member for member in members if member.startswith('RPM[')]
            dna_members = [member for member in members if member.startswith('DPM[')]

            assert len(rna_members) + len(dna_members) == len(members), (
                f'Unexpected cluster members in {cluster_id}: '
                f'{",".join(member for member in members if not member.startswith(("RPM[", "DPM[")))}'
            )

            rna_original = _make_unique(rna_members)
            dna_original = _make_unique(dna_members)

            result.append({
                'cluster_index': cluster_index,
                'cluster_id': cluster_id,
                'n_members': len(members),
                'n_rna_members': len(rna_members),
                'n_dna_members': len(dna_members),
                'rna_original': rna_original,
                'dna_original': dna_original,
                'n_rna_entities': len(rna_original),
                'n_dna_entities': len(dna_original)
            })

    result = pd.DataFrame(result)
    assert (result['n_members'] == result['n_rna_members'] + result['n_dna_members']).all()
    assert result['cluster_index'].is_unique

    return result


def _map_dna_entities2yalid(dna_entities: pd.Index) -> pd.Series:
    result = pd.Series(index=dna_entities, dtype='object')

    if dna_entities.empty:
        return result

    peaks = dna_entities.to_series(index=dna_entities, name='name').to_frame()
    peaks[['chr', 'start', 'end']] = peaks['name'].str.extract(DNA_REGEX)
    assert peaks[['chr', 'start', 'end']].notna().all().all()

    peaks['start'] = peaks['start'].astype('int')
    peaks['end'] = peaks['end'].astype('int')
    peaks['score'] = '.'
    peaks['strand'] = '.'
    peaks = peaks[BED_COLUMNS]

    result = best_left_intersect(
        peaks,
        load_chromhmm_annotation(),
        stranded=False,
        unify_chr_assembly='hg38'
    )
    result = result.set_index('name', verify_integrity=True)['name2']
    result = result.where(result.notna(), None)

    return result


def _map_rna_entities2yagid(rna_entities: pd.Index) -> pd.Series:
    result = pd.Series(index=rna_entities, dtype='object')

    if rna_entities.empty:
        return result

    rna_entities = rna_entities.to_series(index=rna_entities)

    refseq_ids = rna_entities.str.extract(RNA_REFSEQ_REGEX, expand=False)
    has_refseq = refseq_ids.notna()
    if has_refseq.any():
        mapped = id2yagid(refseq_ids[has_refseq])
        mapped = mapped.where(mapped.str.startswith('YAGID'))
        result.loc[mapped.index] = mapped

    unresolved = result.isna()
    if unresolved.any():
        peaks = rna_entities[unresolved].to_frame(name='name')
        peaks[['chr', 'start', 'end']] = peaks['name'].str.extract(RNA_COORD_REGEX)
        peaks = peaks[peaks[['chr', 'start', 'end']].notna().all(axis=1)]

        if not peaks.empty:
            peaks['start'] = peaks['start'].astype('int')
            peaks['end'] = peaks['end'].astype('int')
            peaks['score'] = '.'
            peaks['strand'] = '.'
            peaks = peaks[BED_COLUMNS]

            mapped = best_left_intersect(
                peaks,
                load_extended_annotation(convert2bed=True),
                stranded=False,
                unify_chr_assembly='hg38'
            )
            mapped = mapped.set_index('name', verify_integrity=True)['name2']
            mapped = mapped[mapped.notna()]
            mapped = id2yagid(mapped)
            mapped = mapped.where(mapped.str.startswith('YAGID'))
            result.loc[mapped.index] = mapped

    result = result.where(result.notna(), None)
    return result


def _flatten_unique(column: pd.Series) -> pd.Index:
    values = [item for items in column for item in items]
    return pd.Index(values, dtype='object').unique()


def _apply_entity_mapping(
        original_entities: tuple[str, ...],
        mapping: pd.Series
    ) -> tuple[tuple[str, ...], tuple[str, ...]]:
    mapped_entities = []
    unmapped_entities = []

    for entity in original_entities:
        mapped = mapping.at[entity]
        if mapped is None:
            unmapped_entities.append(entity)
        else:
            mapped_entities.append(mapped)

    return _make_unique(mapped_entities), tuple(unmapped_entities)


@memory.cache
def load_rdsprite_data(filepath_or_url: str = RDSPRITE_DATA_URL) -> pd.DataFrame:
    filepath = _as_local_path(filepath_or_url)
    result = _parse_rdsprite_clusters(filepath)

    dna_mapping = _map_dna_entities2yalid(_flatten_unique(result['dna_original']))
    rna_mapping = _map_rna_entities2yagid(_flatten_unique(result['rna_original']))

    mapped_rna = []
    unmapped_rna = []
    mapped_dna = []
    unmapped_dna = []

    for row in result.itertuples(index=False):
        row_rna, row_rna_unmapped = _apply_entity_mapping(row.rna_original, rna_mapping)
        row_dna, row_dna_unmapped = _apply_entity_mapping(row.dna_original, dna_mapping)

        mapped_rna.append(row_rna)
        unmapped_rna.append(row_rna_unmapped)
        mapped_dna.append(row_dna)
        unmapped_dna.append(row_dna_unmapped)

    result['rna_yagid'] = mapped_rna
    result['rna_unmapped'] = unmapped_rna
    result['dna_yalid'] = mapped_dna
    result['dna_unmapped'] = unmapped_dna

    result['n_yagid'] = result['rna_yagid'].str.len()
    result['n_rna_unmapped'] = result['rna_unmapped'].str.len()
    result['n_yalid'] = result['dna_yalid'].str.len()
    result['n_dna_unmapped'] = result['dna_unmapped'].str.len()

    return result
