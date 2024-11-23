from typing import Callable

import pandas as pd
from tqdm.auto import tqdm


def _load_karr_seq(
        path, *,
        filter_func: Callable = lambda df: df,
        chunksize: int|None = None
    ) -> pd.DataFrame:

    columns = (
        'readID',
        'seqid1', 'pos1',
        'seqid2', 'pos2',
        'strand1', 'strand2'
    )
    kwargs = dict(sep='\t', header=None, names=columns, dtype='str')

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

    assert result['pos1'].str.isdigit().all()
    assert result['pos2'].str.isdigit().all()

    return result


def load_karr_seq(chunksize: int|None, verbose: bool = True) -> pd.DataFrame:
    prefix = 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5064nnn'
    urls = {
        'Nuclear': {
            'R01': f'{prefix}/GSM5064771/suppl/GSM5064771_G1_kethoxal-K562-Nuclear_M17_R01.dedup.pairs.gz',
            'R02': f'{prefix}/GSM5064772/suppl/GSM5064772_G1_kethoxal-K562-Nuclear_M17_R02.dedup.pairs.gz'
        },
        'Total': {
            'R01': f'{prefix}/GSM5064767/suppl/GSM5064767_G1_kethoxal-K562_M15_R01.dedup.pairs.gz',
            'R02': f'{prefix}/GSM5064768/suppl/GSM5064768_G1_kethoxal-K562_M15_R02.dedup.pairs.gz'
        }
    }

    data = []
    for frac in urls:
        for repl in urls[frac]:
            url = urls[frac][repl]
            if verbose:
                print(f'Frac: {frac}; Repl: {repl}')
                print(f'\tDownloading KARR-seq data from {url}')
            data.append(_load_karr_seq(
                url,
                filter_func=lambda df: df[df['seqid1'] != df['seqid2']],
                chunksize=chunksize
            ))
            data[-1]['frac'], data[-1]['repl'] = frac, repl
    return pd.concat(data)
