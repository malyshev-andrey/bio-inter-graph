from urllib.parse import urlencode

import pandas as pd


def _load_encode_metadata(*, cell_line: str|None = None, assay: str, **kwargs) -> pd.DataFrame:
    params = {
        'type': 'File',
        'assay_title': assay,
        'status': 'released'
    }
    if cell_line is not None:
        params['biosample_ontology.term_name'] = cell_line
    params.update(kwargs)

    url = f'https://www.encodeproject.org/report.tsv?{urlencode(params)}'
    metadata = pd.read_csv(url, sep='\t', skiprows=1, dtype='str')

    metadata = metadata.loc[:, ~metadata.isna().all()]

    return metadata


def load_encode_eCLIP(**kwargs) -> pd.DataFrame:
    default_kwargs = dict(
        assay='eCLIP',
        processed='true',
        file_format='bed'
    )
    default_kwargs.update(kwargs)
    metadata = _load_encode_metadata(**default_kwargs)

    assert metadata['Biological replicates'].isin({'1', '2', '1,2'}).all()
    assert metadata['Biological replicates'].value_counts(normalize=True).eq(1/3).all()
    metadata = metadata[metadata['Biological replicates'].eq('1,2')]

    result = {}
    for _, row in tqdm(metadata.iterrows(), total=metadata.shape[0]):
        assembly = row['Genome assembly']
        if assembly not in result:
            result[assembly] = []

        url = f'https://www.encodeproject.org{row["Download URL"]}'
        bed = pd.read_csv(
            url, sep='\t', usecols=range(6), header=None,
            names=['chr', 'start', 'end', 'name', 'score', 'strand']
        )
        bed['name'] = row['Target label']
        bed['cell_line'] = row['Biosample name']
        result[assembly].append(bed)

    for assembly in result:
        result[assembly] = pd.concat(result[assembly])
    return result
