import requests
import pandas as pd
from tqdm.auto import tqdm

from ..shared import memory, CHUNKSIZE


@memory.cache
def _retrieve_karr_seq_metadata(cell_line: str|None = None) -> pd.DataFrame:
    metadata = pd.read_csv(
        'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166155/suppl/filelist.txt',
        sep='\t',
        usecols=['Name']
    )['Name']
    metadata = metadata.set_axis(metadata)
    metadata = metadata[~metadata.eq('GSE166155_RAW.tar')]

    regex = (
        r'^(?P<accession>GSM\d{7})_'
        r'(?P<dendrimers>G\d)_'
        r'(?P<conditions>[^_]+)_'
        r'(?P<group>[BM]\d{2})_'
        r'(?P<repl>R0[12])'
        r'\.dedup\.pairs\.gz$'
    )
    assert metadata.str.match(regex).all()
    metadata = metadata.str.extract(regex)
    assert not metadata.isna().any().any()

    frac_regex = r'-(?P<frac>Total|Nuclear)(RNA)?$'
    metadata['frac'] = metadata['conditions'].str.extract(frac_regex)['frac']
    metadata['frac'] = metadata['frac'].fillna('Total')

    cell_line_regex = '^kethoxal-(?P<cell_line>[^-+]+)(\+S2)?(-.*)?$'
    metadata['cell_line'] = metadata['conditions'].str.extract(cell_line_regex)['cell_line']

    in_vivo_conditions = [
        'kethoxal-F123',
        'kethoxal-HepG2-TotalRNA',
        'kethoxal-K562-Nuclear',
        'kethoxal-mESC',
        'kethoxal-K562',
        'kethoxal-HepG2',
        'kethoxal-HEK293T'
    ]

    metadata['is_in_vivo'] = metadata['conditions'].isin(in_vivo_conditions)
    metadata = metadata[~metadata['cell_line'].isna()]
    metadata = metadata[~metadata['cell_line'].eq('riboplus')]

    if cell_line is not None:
        cell_lines = metadata['cell_line'].unique()
        if cell_line not in cell_lines:
            raise ValueError(
                f'"{cell_line}" is not a valid argument. '
                f'Valid arguments are: {", ".join(cell_lines)}'
            )
        metadata = metadata[metadata['cell_line'].eq(cell_line)]

    metadata = metadata[metadata['is_in_vivo']]

    metadata['url'] = metadata.apply(
        lambda row: (
            f'https://ftp.ncbi.nlm.nih.gov/geo/samples'
            f'/{row["accession"][:-3] + "nnn"}/{row["accession"]}/suppl/{row.name}'
        ),
        axis=1
    )
    for url in metadata['url']:
        requests.head(url, allow_redirects=True, timeout=5).raise_for_status()

    metadata['url'] = metadata['url'].str.replace(r'https://', 'ftp://')
    return metadata


def _load_single_karr_seq(
        path, *,
        filter_func: Callable = lambda df: df,
        chunksize: int|None = CHUNKSIZE
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

    assert 'pos1' not in result.columns or result['pos1'].str.isdigit().all()
    assert 'pos2' not in result.columns or result['pos2'].str.isdigit().all()

    return result
