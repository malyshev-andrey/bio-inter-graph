import pandas as pd
from tqdm.auto import tqdm

GFF_COLUMNS = [
    'chr', 'source', 'type',
    'start', 'end', 'score',
    'strand', 'phase', 'attributes'
]

def read_feature_table(
        path: str, *,
        filter_func=lambda df: df,
        chunksize: int|None = None
    ) -> pd.DataFrame:

    kwargs = dict(
        sep='\t', comment='#',
        header=None, names=GFF_COLUMNS,
        dtype='str'
    )

    if chunksize is None:
        return filter_func(pd.read_csv(path, **kwargs))
    else:
        result = []
        desc = path if len(path) < 40 else path[:20] + ' ... ' + path[-20:]
        with tqdm(desc=desc) as progress_bar:
            for chunk in pd.read_csv(path, chunksize=chunksize, **kwargs):
                progress_bar.update(chunk.shape[0])
                result.append(filter_func(chunk))
        return pd.concat(result)
