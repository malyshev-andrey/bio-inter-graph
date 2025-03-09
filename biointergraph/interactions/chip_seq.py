from urllib.parse import urlencode

import pandas as pd
from tqdm.auto import tqdm

from ..annotations import bed_merge, bed_intersect, sanitize_bed, load_ChromHMM_annotation
from ..ids_mapping import id2yapid
from ..shared import memory, BED_COLUMNS


def _load_encode_ChIP_seq_metadata(cell_line: str|None = None) -> pd.DataFrame:

    params = [
        ('assay_title', 'TF ChIP-seq'),
        ('status', 'released'),
        ('file_format', 'bed'),
        ('assembly', 'GRCh38'),
        ('output_type', 'IDR thresholded peaks'),
        ('output_type', 'conservative IDR thresholded peaks')
    ]
    if cell_line is not None:
        params.append(('biosample_ontology.term_name', cell_line))
    params = urlencode(params)

    url = f'https://www.encodeproject.org/report.tsv?type=File&{params}'
    metadata = pd.read_csv(url, sep='\t', skiprows=1, dtype='str')

    n_experiments = metadata['Dataset'].nunique()

    metadata['is_conservative'] = metadata['Output type'].eq('conservative IDR thresholded peaks')
    metadata = metadata[
        ~metadata.groupby(['Dataset'])['is_conservative'].transform('any')
        | metadata['is_conservative']
    ]

    assert metadata['Dataset'].nunique() == n_experiments

    metadata['Date created'] = pd.to_datetime(metadata['Date created'])
    latest_file_date = metadata.groupby('Dataset')['Date created'].transform('max')
    metadata = metadata[metadata['Date created'] == latest_file_date]
    assert metadata['Dataset'].is_unique

    return metadata


@memory.cache
def load_encode_chip_seq_peaks(cell_line: str|None = None) -> pd.DataFrame:
    metadata = _load_encode_ChIP_seq_metadata(cell_line)

    result = []
    with tqdm(desc='Peaks:') as progress_bar:
        for _, row in tqdm(metadata.iterrows(), total=metadata.shape[0], desc='Experiments:'):
            bed = pd.read_csv(
                f'https://www.encodeproject.org{row["Download URL"]}',
                sep='\t', usecols=range(6),
                header=None, names=BED_COLUMNS,
                dtype='str'
            )
            bed['name'] = row['Target label']
            result.append(bed)
            progress_bar.update(bed.shape[0])
    result = pd.concat(result)
    result = sanitize_bed(result, stranded=False)

    result = bed_merge(result, by='Name')
    result['score'], result['strand'] = 1000, '.'

    return result


@memory.cache
def load_chip_seq_data() -> pd.DataFrame:
    peaks = load_encode_chip_seq_peaks('K562')
    ChromHMM = load_ChromHMM_annotation()
    peaks['name'] = id2yapid('SYMBOL:' + peaks['name'])
    mapped = peaks['name'].str.startswith('YAPID')
    unmapped = peaks[~mapped]['name'].unique()
    print('ChIP-seq data unmapped protein symbols:', *unmapped)
    print(f'Invalid peaks frac: {1-mapped.mean():.04f}')

    peaks = peaks[mapped]
    result = bed_intersect(ChromHMM, peaks, strandedness=None, unify_chr_assembly='hg38', jaccard=True)
    peak_id = ['chr', 'start2', 'end2', 'name2']
    result = result.sort_values('jaccard', ascending=False)
    result = result.drop_duplicates(peak_id, keep='first')
    return result
