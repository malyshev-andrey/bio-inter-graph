from io import StringIO
from math import ceil
from time import sleep
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from ..shared import GFF_COLUMNS
from ..ids import drop_id_version


def refseq_transcript_id2entrez_gene_id(ids: pd.Series, chunksize: int = 10000) -> pd.Series:
    unique_ids = ids.unique()
    n_ids = unique_ids.shape[0]
    n_chunks = ceil(n_ids / chunksize)
    params = dict(
        db="nuccore",
        retmode='text',
        rettype='gff3'
    )
    ids_data = []
    with tqdm(desc='refseq2entrez IDs processed:', total=n_ids) as progress_bar:
        with ThreadPoolExecutor(max_workers=n_chunks) as executor:
            futures = []
            for chunk in np.array_split(unique_ids, n_chunks):
                params['id'] = ','.join(chunk)
                futures.append(executor.submit(
                    requests.post,
                    'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi',
                    data=params
                ))
                sleep(1)
            for future in as_completed(futures):
                response = future.result()
                response.raise_for_status()
                data = pd.read_csv(
                    StringIO(response.text),
                    comment='#',
                    sep='\t',
                    header=None,
                    names=GFF_COLUMNS
                )
                progress_bar.update(data['chr'].nunique())
                ids_data.append(data)

    ids_data = pd.concat(ids_data)

    ids_data['chr'] = drop_id_version(ids_data['chr'])

    regex = r'GeneID:(?P<gene_id>\d+)'
    assert (ids_data['attributes'].str.count(regex) <= 1).all()
    ids_data['gene_id'] = ids_data['attributes'].str.extract(regex)['gene_id']
    ids_data = ids_data[~ids_data['gene_id'].isna()]
    ids_data = ids_data[['chr', 'gene_id']].drop_duplicates()

    assert unique_ids.shape[0] == ids_data.shape[0]

    ids_data = ids_data['gene_id'].set_axis(ids_data['chr'])
    result = ids.map(ids_data)
    return result
