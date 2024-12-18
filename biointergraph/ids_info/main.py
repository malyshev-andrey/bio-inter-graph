from io import StringIO
from math import ceil

import requests
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from ..shared import GFF_COLUMNS


def retrieve_gene_id4refseq_transcripts(ids: pd.Series, *, chunksize: int = 250, **kwargs) -> pd.Series:
    unique_ids = ids.unique()
    n_ids = unique_ids.shape[0]
    n_chunks = ceil(n_ids / chunksize)

    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'

    params = dict(
        db="nuccore",
        retmode='text',
        rettype='gff3',
        **kwargs
    )

    ids_data = []

    with tqdm(desc='IDs processed:', total=n_ids) as progress_bar:
        for chunk in np.array_split(unique_ids, n_chunks):
            progress_bar.update(chunk.shape[0])

            params['id'] = ','.join(chunk)
            resp = requests.post(url, params=params)
            resp.raise_for_status()
            resp = StringIO(resp.text)

            ids_data.append(pd.read_csv(
                resp, sep='\t', comment='#', header=None, names=GFF_COLUMNS
            ))

    ids_data = pd.concat(ids_data)

    ids_data['chr'] = ids_data['chr'].str.split('.', expand=True)[0]

    regex = r'GeneID:(\d+)'
    assert (ids_data['attributes'].str.count(regex) <= 1).all()
    ids_data['gene_id'] = ids_data['attributes'].str.extract(regex)[0]
    ids_data = ids_data[~ids_data['gene_id'].isna()]
    ids_data = ids_data[['chr', 'gene_id']].drop_duplicates()

    assert unique_ids.shape[0] == ids_data.shape[0]

    ids_data = ids_data.set_index('chr', verify_integrity=True)

    ids_data = ids_data.loc[ids, 'gene_id']
    ids_data.index = ids.index

    return ids_data
