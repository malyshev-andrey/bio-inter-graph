from typing import Callable

import requests
from requests.exceptions import HTTPError
import pandas as pd
from tqdm.auto import tqdm

from ..ids import drop_id_version


def _retrieve_ucsc_schema(table, assembly: str = 'hg38') -> list[str]:
    assert assembly in ['hg19', 'hg38']
    url = f'https://api.genome.ucsc.edu/list/schema?genome={assembly};track={table}'
    response = requests.get(url)
    try:
        response.raise_for_status()
    except HTTPError:
        if table == 'chromAlias':
            return ['alias', 'chrom', 'source']
        else:
            raise
    response = response.json()
    assert 'columnTypes' in response, f'Failed to retrieve schema from {url}'
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


def unify_chr(chr: pd.Series, assembly: str = 'hg38') -> pd.Series:
    """
    Unify chromosome names using UCSC chromAlias table.

    This function standardizes chromosome names in a pandas Series using
    the UCSC `chromAlias` table for the specified genome assembly. It maps
    identifiers to their canonical UCSC names and attempts to handle identifiers
    without version suffixes. If no mapping is found, the original name is
    retained in the output.

    Args:
        chr (pd.Series): A pandas Series containing chromosome names to be unified.
        assembly (str): Genome assembly name. Supported assemblies include:
            - "GRCh38" or "hg38" (default)
            - "GRCh37" or "hg19"

    Returns:
        pd.Series: A pandas Series with unified chromosome names.

    Raises:
        ValueError: If an invalid assembly name is provided.

    Notes:
        - The function downloads the UCSC chromAlias table directly from
            UCSC servers
        - The function uses the `ids.drop_id_version` function to strip version
            suffixes from chromosome identifiers for better mapping.
    """

    ASSEMBLIES = {
        'GRCh38': 'hg38',
        'hg38': 'hg38',
        'GRCh37': 'hg19',
        'hg19': 'hg19'
    }
    if assembly not in ASSEMBLIES:
        raise ValueError(
            f'"{assembly}" is not a valid argument. '
            f'Valid arguments are: {", ".join(ASSEMBLIES)}'
        )
    assembly = ASSEMBLIES[assembly]

    mapping = fetch_ucsc_table('chromAlias', assembly=assembly)
    mapping = mapping.set_index('alias', verify_integrity=True)['chrom']
    drop_version_map = drop_id_version(chr).map(mapping)

    result = chr.map(mapping).combine_first(drop_version_map).combine_first(chr)

    return result
