import os
import stat
import tempfile
import subprocess
import shlex
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import requests
from tqdm.auto import tqdm

from ..shared import BED_COLUMNS, _read_tsv, memory, remote_file2local
from .main import _annotate_peaks
from ..annotations import load_chromhmm_annotation, sanitize_bed
from ..ids_mapping import id2yapid


def _bigbed2bed(path_or_url: str, name: str, *, converter: str) -> pd.DataFrame:
    bed = tempfile.NamedTemporaryFile(delete=False)
    bed.close()

    path_or_url = remote_file2local(path_or_url, progress_bar=False)

    cmd = f'{converter} {path_or_url} {bed.name}'
    subprocess.run(shlex.split(cmd), check=True)

    result = _read_tsv(
        bed.name,
        header=None,
        usecols=range(6),
        names=BED_COLUMNS,
        chunksize=None
    )
    result['name'] = name

    os.remove(bed.name)

    return result


def _gtrd_metadata2bed(metadata: pd.DataFrame) -> pd.DataFrame:
    converter = tempfile.NamedTemporaryFile(delete=False)
    os.chmod(
        converter.name,
        os.stat(converter.name).st_mode | stat.S_IEXEC
    )
    response = requests.get('https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed')
    response.raise_for_status()
    converter.write(response.content)
    converter.flush()
    converter.close()

    result = []
    with ThreadPoolExecutor(max_workers=100) as executor:
        futures = []
        for _, row in metadata.iterrows():
            futures.append(executor.submit(
                _bigbed2bed,
                f'http://gtrd.biouml.org:8888{row["path"]}',
                row["uniprot"],
                converter=converter.name
            ))

        tqdm_kwargs = dict(total=len(futures), unit='file', desc='GTRD ChIP-seq')
        for future in tqdm(as_completed(futures), **tqdm_kwargs):
            result.append(future.result())
    result = pd.concat(result)

    os.remove(converter.name)

    result['score'], result['strand'] = 1000, '.'
    result = sanitize_bed(result, stranded=False)
    return result


def _load_gtrd_metadata(cell_line: str|None = None) -> pd.DataFrame:
    result = pd.read_csv(
        'http://gtrd.biouml.org:8888/downloads/current/bigBeds/hg38/ChIP-seq/Meta-clusters_by_TF_and_Cell_Type',
        sep='\t',
        header=None,
        skiprows=14,
        skipfooter=4,
        names=['html'],
        engine='python'
    )
    result['path'] = result['html'].str.extract(r"href='([^']+)'", expand=False)
    result['path'] = result['path'].str.replace('/egrid', '/downloads/current')
    result['file'] = result['path'].str.split('/', expand=True).iloc[:,-1]

    regex = r'^(?P<symbol>[A-Z0-9]+)_(?P<uniprot>[A-Z0-9]{6})_Meta-clusters_(?P<cell_id>\d+).bb$'
    assert result['file'].str.match(regex).all()
    result = pd.concat([
        result,
        result['file'].str.extract(regex)
    ], axis=1)

    cell_types = _read_tsv(
        'http://gtrd.biouml.org:8888/downloads/current/metadata/cell_types_and_tissues.metadata.txt',
        chunksize=None,
        usecols=['id', 'title', 'species']
    )
    cell_types = cell_types.rename(columns={'id': 'cell_id'})

    cell_types['cell_line'] = cell_types['title'].str.split(' ', expand=True)[0]

    result = result.merge(cell_types, how='left', validate='many_to_one')

    if cell_line is not None:
        result = result[result['cell_line'].eq(cell_line)]

    return result


@memory.cache
def load_gtrd_chip_seq_data(cell_line: str|None = None) -> pd.DataFrame:
    metadata = _load_gtrd_metadata(cell_line)
    peaks = _gtrd_metadata2bed(metadata)

    annotation = load_chromhmm_annotation()

    result = _annotate_peaks(
        peaks, annotation,
        assembly='hg38',
        stranded=False,
        convert_ids=False
    )

    result['source'] = id2yapid(result['source'], strict=True)
    result = result.drop_duplicates()

    return result
