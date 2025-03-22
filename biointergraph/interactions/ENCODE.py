from urllib.parse import urlencode
from typing import Iterable

import pandas as pd
from tqdm.auto import tqdm

from ..shared import memory, BED_COLUMNS
from ..annotations import load_refseq_bed, load_gencode_bed, sanitize_bed, bed_intersect
from ..ids_mapping import id2yapid, id2yagid


def load_encode_metadata(
        assay: str|Iterable[str] = (), *,
        entity_type: str = 'File',
        cell_line: str|None = None,
        released: bool = True,
        **kwargs
    ) -> pd.DataFrame:

    params = []
    if isinstance(assay, str):
        assay = [assay]
    params.extend(('assay_title', title) for title in assay)

    if cell_line is not None:
        params.append(('biosample_ontology.term_name', cell_line))

    if released:
        params.append(('status', 'released'))

    params.extend(kwargs.items())
    params = urlencode(params)

    url = f'https://www.encodeproject.org/report.tsv?type={entity_type}&{params}'
    print(f'ENCODE metadata URL: {url}')
    metadata = pd.read_csv(url, sep='\t', skiprows=1, dtype='str')

    metadata = metadata.loc[:, ~metadata.isna().all()]
    assert metadata.shape[0] > 0, 'No metadata found!'

    return metadata


@memory.cache
def _load_encode_eCLIP(assembly: str, cell_line: str|None = None) -> pd.DataFrame:
    ASSEMBLIES = {
        'hg38': 'GRCh38', 'GRCh38': 'GRCh38',
        'GRCh37': 'hg19', 'hg19': 'hg19',
    }
    if assembly not in ASSEMBLIES:
        raise ValueError(
            f'"{assembly}" is not a valid argument. '
            f'Valid arguments are: {", ".join(ASSEMBLIES)}'
        )
    assembly = ASSEMBLIES[assembly]

    default_kwargs = dict(
        assay='eCLIP',
        processed='true',
        file_format='bed',
        assembly=assembly
    )
    if cell_line is not None:
        default_kwargs['cell_line'] = cell_line
    metadata = load_encode_metadata(**default_kwargs)

    replicates = metadata['Biological replicates']
    assert replicates.isin({'1', '2', '1,2'}).all()
    assert replicates.value_counts(normalize=True).eq(1/3).all()
    metadata = metadata[replicates.eq('1,2')]

    result = []
    with tqdm(desc='ENCODE eCLIP', unit='peak') as progress_bar:
        for _, row in tqdm(metadata.iterrows(), total=metadata.shape[0], desc='ENCODE eCLIP', unit='experiment'):
            bed = pd.read_csv(
                f'https://www.encodeproject.org{row["Download URL"]}',
                sep='\t', usecols=range(6),
                header=None, names=BED_COLUMNS,
                dtype='str'
            )
            bed['name'] = row['Target label']
            bed['cell_line'] = row['Biosample name']

            result.append(bed)
            progress_bar.update(bed.shape[0])

    result = pd.concat(result)
    result = sanitize_bed(result)

    return result


@memory.cache
def encode_eCLIP2pairwise(
        assembly: str,
        annotation: str,
        cell_line: str|None = None
    ) -> pd.DataFrame:
    eCLIP_bed = _load_encode_eCLIP(assembly=assembly, cell_line=cell_line)
    annotation_bed = {
        'gencode': load_gencode_bed,
        'refseq': load_refseq_bed
    }[annotation](assembly=assembly, feature='gene')

    result = bed_intersect(
        eCLIP_bed,
        annotation_bed,
        unify_chr_assembly=assembly,
        jaccard=True,
        how='left'
    )

    no_intersect = result['start2'].eq(-1)
    print(f'ENCODE eCLIP peaks without intersections: {no_intersect.sum()}')
    result = result[~no_intersect]

    peak_id = {f'{c}1': c for c in BED_COLUMNS}
    result = result.rename(columns=peak_id)
    result = result.sort_values('jaccard')
    result = result.drop_duplicates(peak_id.values(), keep='last')

    result['yapid'] = id2yapid('SYMBOL:' + result['name'])
    assert result['yapid'].str.startswith('YAPID').all()
    result['yagid'] = id2yagid(result['name2'])
    assert result['yagid'].str.startswith('YAGID').all()

    result = result[['yapid', 'yagid']]
    result = result.drop_duplicates()

    return result
