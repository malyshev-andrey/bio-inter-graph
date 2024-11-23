from warnings import warn

import pandas as pd
from tqdm.auto import tqdm

from ..shared import GFF_COLUMNS


def _validate_feature_table(ft: pd.DataFrame) -> pd.DataFrame:
    is_valid_start = ft['start'].str.isdigit()
    is_valid_end = ft['end'].str.isdigit()
    is_valid_strand = ft['strand'].isin({'+', '-', '.'})
    if not is_valid_start.all():
        incorrect = " ".join(ft['start'][~is_valid_start].head(5))
        warn(
            f'Some start coordinates seem to be incorrect: {incorrect}.\n'
            'Skipping these lines!'
        )
    if not is_valid_end.all():
        incorrect = " ".join(ft['end'][~is_valid_end].head(5))
        warn(
            f'Some end coordinates seem to be incorrect: {incorrect}.\n'
            'Skipping these lines!'
        )
    if not is_valid_strand.all():
        incorrect = " ".join(ft['strand'][~is_valid_strand].head(5))
        warn(
            f'Some strand values seem to be incorrect: {incorrect}.\n'
            'Skipping these lines!'
        )
    ft = ft[is_valid_strand & is_valid_end & is_valid_start]
    return ft


def read_feature_table(
        path: str, *,
        filter_func=lambda df: df,
        chunksize: int|None = None,
        validation: bool = True
    ) -> pd.DataFrame:

    kwargs = dict(
        sep='\t', comment='#',
        header=None, names=GFF_COLUMNS,
        dtype='str'
    )

    if chunksize is None:
        result = filter_func(pd.read_csv(path, **kwargs))
    else:
        result = []
        desc = path if len(path) < 40 else path[:20] + ' ... ' + path[-20:]
        with tqdm(desc=desc) as progress_bar:
            for chunk in pd.read_csv(path, chunksize=chunksize, **kwargs):
                progress_bar.update(chunk.shape[0])
                result.append(filter_func(chunk))
        result = pd.concat(result)

    if validation:
        result = _validate_feature_table(result)

    result['start'] = result['start'].astype('int')
    result['end'] = result['end'].astype('int')

    return result
