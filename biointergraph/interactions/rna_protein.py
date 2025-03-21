import pandas as pd

from ..shared import _read_tsv, memory


def _load_postar3_peaks(species: str = 'human', **kwargs) -> pd.DataFrame:
    default_kwargs = dict(
        header=None,
        names=[
            'chr', 'start', 'end', 'peak_id', 'strand', 'name',
            'method', 'cell_line', 'accessions', 'score'
        ]
    )
    default_kwargs.update(kwargs)

    result = _read_tsv(
        f'https://cloud.tsinghua.edu.cn/d/8133e49661e24ef7a915/files/?dl=1&p={species}.txt.gz',
        **default_kwargs
    )
    return result


@memory.cache
def load_postar3_data(species: str, cell_line: str, **kwargs) -> pd.DataFrame:
    result = _load_postar3_peaks(
        species=species,
        filter_func=lambda df: df[df['cell_line'].eq(cell_line)],
        **kwargs
    )
    return result
