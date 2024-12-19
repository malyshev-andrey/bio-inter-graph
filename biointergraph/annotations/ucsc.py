from typing import Callable

import requests
import pandas as pd
from tqdm.auto import tqdm


def _retrieve_ucsc_schema(table, assembly: str = 'hg38') -> list[str]:
    response = requests.get(
        f'https://api.genome.ucsc.edu/list/schema?genome={assembly};track={table}'
    )
    response = response.json()
    response = response['columnTypes']
    response = [column['name'] for column in response]

    return response


def fetch_ucsc_table(
        table,
        assembly: str = 'hg38',
        chunksize: int|None = None,
        filter_func: Callable[[pd.DataFrame], pd.DataFrame] = lambda df: df,
        **kwargs
    ) -> pd.DataFrame:
    default_kwargs = dict(
        names=_retrieve_ucsc_schema(table, assembly),
        header=None,
        sep='\t'
    )
    default_kwargs.update(kwargs)

    url = f'https://hgdownload.soe.ucsc.edu/goldenPath/{assembly}/database/{table}.txt.gz'

    if chunksize is None:
        result = filter_func(pd.read_csv(url, **default_kwargs))
    else:
        result = []
        with tqdm(desc=url) as progress_bar:
            for chunk in pd.read_csv(url, chunksize=chunksize, **default_kwargs):
                progress_bar.update(chunk.shape[0])
                result.append(filter_func(chunk))
        result = pd.concat(result)
    return result
