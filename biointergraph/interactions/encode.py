from urllib.parse import urlencode
from typing import Iterable, Callable
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from tqdm.auto import tqdm

from ..shared import memory, BED_COLUMNS, _read_tsv
from ..annotations import (
    load_refseq_bed, load_gencode_bed,
    sanitize_bed, bed_cluster,
    load_chromhmm_annotation
)
from ..ids_mapping import id2yapid
from .main import _annotate_peaks


IDR_BED_COLUMNS = BED_COLUMNS + [
    'signalValue',
    'p_value', 'q_value',
    'summit',
    'localIDR', 'globalIDR'
]

def load_encode_metadata(
        assay: str|Iterable[str] = (), *,
        entity_type: str = 'File',
        cell_line: str|None = None,
        released: bool = True,
        **kwargs
    ) -> pd.DataFrame:

    params = []
    if isinstance(assay, str):
        assay = [assay]
    params.extend(('assay_title', title) for title in assay)

    if cell_line is not None:
        params.append(('biosample_ontology.term_name', cell_line))

    if released:
        params.append(('status', 'released'))

    if 'assembly' in kwargs:
        assembly = kwargs['assembly']
        ASSEMBLIES = {
            'hg38': 'GRCh38', 'GRCh38': 'GRCh38',
            'GRCh37': 'hg19', 'hg19': 'hg19',
        }
        if assembly not in ASSEMBLIES:
            raise ValueError(
                f'"{assembly}" is not a valid argument. '
                f'Valid arguments are: {", ".join(ASSEMBLIES)}'
            )
        kwargs['assembly'] = ASSEMBLIES[assembly]

    params.extend(kwargs.items())
    params = urlencode(params)

    url = f'https://www.encodeproject.org/report.tsv?type={entity_type}&{params}'
    print(f'ENCODE metadata URL: {url}')
    metadata = pd.read_csv(url, sep='\t', skiprows=1, dtype='str')

    metadata = metadata.loc[:, ~metadata.isna().all()]
    assert metadata.shape[0] > 0, 'No metadata found!'

    return metadata


def _encode_metadata_row2bed(
        row: pd.Series,
        features: str|dict|Iterable[str]|None = None,
        colnames: list[str] = BED_COLUMNS,
        filter_func: Callable[[pd.DataFrame], pd.DataFrame] = lambda df: df
    ) -> pd.DataFrame:
    bed = _read_tsv(
        f'https://www.encodeproject.org{row["Download URL"]}',
        usecols=range(len(colnames)),
        header=None,
        names=colnames,
        chunksize=None,
        filter_func=filter_func
    )
    bed['name'] = row['Target label']

    if isinstance(features, dict):
        for key, value in features.items():
            bed[key] = row[value]
    elif features is not None:
        for name in features:
            bed[name] = row[name]

    return bed


def _encode_metadata2bed(
        files: pd.DataFrame, *,
        features: str|dict|Iterable[str]|None = None,
        desc: str|None = None,
        stranded: bool = True,
        colnames: list[str] = BED_COLUMNS,
        filter_func: Callable[[pd.DataFrame], pd.DataFrame] = lambda df: df
    ) -> pd.DataFrame:
    if desc is None:
        assay = files['Assay term name'].unique().item()
        desc = f'ENCODE {assay}'


    result = []
    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = []
        for _, row in files.iterrows():
            futures.append(executor.submit(
                _encode_metadata_row2bed,
                row, features, colnames, filter_func
            ))

        tqdm_kwargs = dict(desc=desc, total=len(futures), unit='file')
        for future in tqdm(as_completed(futures), **tqdm_kwargs):
            result.append(future.result())

    result = pd.concat(result)

    result = sanitize_bed(result, stranded=stranded)
    return result


@memory.cache
def _load_encode_eclip_bed(assembly: str, cell_line: str|None = None) -> pd.DataFrame:
    default_kwargs = dict(
        assay='eCLIP',
        processed='true',
        file_format='bed',
        assembly=assembly
    )
    if cell_line is not None:
        default_kwargs['cell_line'] = cell_line
    metadata = load_encode_metadata(**default_kwargs)

    replicates = metadata['Biological replicates']
    assert replicates.isin({'1', '2', '1,2'}).all()
    assert replicates.value_counts(normalize=True).eq(1/3).all()
    metadata = metadata[replicates.eq('1,2')]

    result = _encode_metadata2bed(metadata, features={'cell_line': 'Biosample name'})

    return result


@memory.cache
def load_encode_eclip_data(
        assembly: str,
        annotation: str,
        cell_line: str|None = None
    ) -> pd.DataFrame:
    peaks = _load_encode_eclip_bed(assembly=assembly, cell_line=cell_line)
    annotation = {
        'gencode': load_gencode_bed,
        'refseq': load_refseq_bed
    }[annotation](assembly=assembly, feature='gene')

    result = _annotate_peaks(
        peaks, annotation,
        assembly=assembly,
        desc='ENCODE eCLIP',
        convert_ids=True
    )

    return result


@memory.cache
def load_encode_iclip_data(annotation: str, *, cell_line: str|None = None) -> pd.DataFrame:
    metadata = load_encode_metadata(
        assay='iCLIP',
        cell_line=cell_line,
        file_format='bed',
        processed='true',
        assembly='hg19'
    )
    peaks = _encode_metadata2bed(metadata, features={'repl': 'Biological replicates'})

    annotation = {
        'gencode': load_gencode_bed,
        'refseq': load_refseq_bed
    }[annotation](assembly='hg19', feature='gene')

    assert peaks['repl'].nunique() == 2

    result = []
    for _, repl in peaks.groupby('repl'):
        result.append(
            _annotate_peaks(repl, annotation, assembly='hg19', convert_ids=True)
        )
    result = result[0].merge(result[1], how='inner', validate='one_to_one')

    return result


@memory.cache
def load_encode_rip_data(annotation: str, *, cell_line: str|None = None):
    metadata = load_encode_metadata(
        ['RIP-seq', 'RIP-chip'],
        cell_line=cell_line,
        file_format='bed',
        processed='true',
        assembly='hg19'
    )
    metadata = metadata[
        ~metadata['Target label'].isna() &
        ~metadata['Target label'].eq('T7')
    ]

    result = _encode_metadata2bed(metadata, stranded=False, desc='RIP-seq, RIP-chip')

    annotation = {
        'gencode': load_gencode_bed,
        'refseq': load_refseq_bed
    }[annotation](assembly='hg19', feature='gene')

    result = _annotate_peaks(result, annotation, assembly='hg19', stranded=False, convert_ids=True)

    return result


@memory.cache
def _load_encode_chip_seq_bed(assembly: str, cell_line: str|None = None) -> pd.DataFrame:
    metadata = load_encode_metadata(
        'TF ChIP-seq',
        cell_line=cell_line,
        processed='true',
        file_format='bed',
        assembly=assembly
    )

    metadata = metadata[
        metadata['Output type'].eq('IDR ranked peaks') &
        metadata['Analysis step version'].eq('/analysis-step-versions/tf-chip-seq-replicated-idr-step-v-1-0/')
    ]
    metadata = metadata.sort_values('File size')
    metadata = metadata.drop_duplicates('Dataset', keep='last')

    metadata = metadata[
        ~metadata['Target label'].isna() &
        ~metadata['Target label'].isin({'POLR2AphosphoS5', 'POLR2AphosphoS2'})
    ]

    result = _encode_metadata2bed(
        metadata,
        stranded=False,
        colnames=IDR_BED_COLUMNS,
        filter_func=lambda df: df[10**(-df['globalIDR'].astype('float')) < 0.05]
    )

    result = bed_cluster(result, by='Name')

    result['globalIDR'] = result['globalIDR'].astype('float')
    result = result.groupby(['chr', 'name', 'Cluster'], as_index=False, observed=True).agg(
        start=('start', 'min'),
        end=('end', 'max'),
        weight=('globalIDR', 'max')
    )

    result = result.drop(columns='Cluster')

    return result


@memory.cache
def load_encode_chip_seq_data(assembly: str, cell_line: str|None = None) -> pd.DataFrame:
    peaks = _load_encode_chip_seq_bed(assembly, cell_line=cell_line)
    annotation = load_chromhmm_annotation()

    result = _annotate_peaks(
        peaks, annotation,
        assembly='hg38',
        desc='ENCODE ChIP-seq',
        stranded=False,
        drop_duplicates=False
    )

    result['source'] = id2yapid('SYMBOL:' + result['source'], strict=True)

    result = result.groupby(['source', 'target'], as_index=False, observed=True)['weight'].max()

    return result
