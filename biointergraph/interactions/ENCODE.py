from urllib.parse import urlencode

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from .main import summarize_pairwise
from ..shared import BED_COLUMNS
from ..annotations import load_refseq_bed, load_gencode_bed, sanitize_bed, bed_intersect
from ..ids import drop_id_version


def load_encode_metadata(*, cell_line: str|None = None, assay: str, **kwargs) -> pd.DataFrame:
    """
    Loads metadata for ENCODE project files based on the specified cell line and assay.

    This function queries the ENCODE database for metadata related to the specified assay (experimental protocol)
    and optionally filters by a given cell line. It returns the result as a pandas DataFrame containing the metadata
    for released files that match the query parameters.

    Parameters:
    -----------
    cell_line : str or None, optional
        The cell line to filter the metadata by. If None, no filtering by cell line is applied.

    assay : str
        The assay (experimental protocol) to filter the metadata by. This is a required parameter.

    **kwargs : keyword arguments
        Additional parameters to pass as filters for the ENCODE metadata query. These will be appended to the query
        parameters in the URL (e.g., for additional filtering based on specific file types or statuses).

    Returns:
    --------
    pd.DataFrame
        A pandas DataFrame containing the metadata of files retrieved from the ENCODE project database, with columns
        for file attributes such as file format, file size, and more. The dataframe will not contain fully empty columns.

    Raises:
    -------
    AssertionError
        If no files are found matching the query parameters, an assertion error is raised.

    Example:
    --------
    >>> load_encode_metadata(assay="eCLIP", cell_line="K562")
    This will return metadata for the eCLIP assays conducted on the K562 cell line.
    """

    params = {
        'type': 'File',
        'assay_title': assay,
        'status': 'released'
    }
    if cell_line is not None:
        params['biosample_ontology.term_name'] = cell_line
    params.update(kwargs)

    url = f'https://www.encodeproject.org/report.tsv?{urlencode(params)}'
    print(f'ENCODE metadata URL: {url}')
    metadata = pd.read_csv(url, sep='\t', skiprows=1, dtype='str')

    metadata = metadata.loc[:, ~metadata.isna().all()]

    assert metadata.shape[0] > 0, 'No files found!'

    return metadata


def load_encode_eCLIP(assembly: str, cell_line: str|None = None, **kwargs) -> pd.DataFrame:
    ASSEMBLIES = {
        'hg38': 'GRCh38', 'GRCh38': 'GRCh38',
        'GRCh37': 'hg19', 'hg19': 'hg19',
    }
    if assembly not in ASSEMBLIES:
        raise ValueError(
            f'"{assembly}" is not a valid argument. '
            f'Valid arguments are: {", ".join(ASSEMBLIES)}'
        )
    assembly = ASSEMBLIES[assembly]

    default_kwargs = dict(
        assay='eCLIP',
        processed='true',
        file_format='bed',
        assembly=assembly
    )
    if cell_line is not None:
        default_kwargs['cell_line'] = cell_line
    default_kwargs.update(kwargs)
    metadata = load_encode_metadata(**default_kwargs)

    replicates = metadata['Biological replicates']
    assert replicates.isin({'1', '2', '1,2'}).all()
    assert replicates.value_counts(normalize=True).eq(1/3).all()
    metadata = metadata[replicates.eq('1,2')]

    result = []
    with tqdm(desc='eCLIP peaks found') as progress_bar:
        for _, row in tqdm(metadata.iterrows(), total=metadata.shape[0], desc='eCLIP experiments'):
            bed = pd.read_csv(
                f'https://www.encodeproject.org{row["Download URL"]}',
                sep='\t', usecols=range(6),
                header=None, names=BED_COLUMNS,
                dtype='str'
            )
            bed['name'] = row['Target label']
            bed['cell_line'] = row['Biosample name']

            result.append(bed)
            progress_bar.update(bed.shape[0])

    result = pd.concat(result)
    result = sanitize_bed(result)

    return result


def encode_eCLIP2pairwise(assembly: str, annotation: str, cell_line: str|None = None, **kwargs) -> pd.DataFrame:
    eCLIP_bed = load_encode_eCLIP(assembly=assembly, cell_line=cell_line)
    annotation_bed = {
        'gencode': load_gencode_bed,
        'refseq': load_refseq_bed
    }[annotation](assembly=assembly, feature='gene')

    result = bed_intersect(
        eCLIP_bed,
        annotation_bed,
        unify_chr_assembly=assembly,
        **kwargs
    )
    intersect = (
        np.minimum(result['end1'], result['end2'])
        - np.maximum(result['start1'], result['start2'])
    )
    covered_peak_frac = intersect / (result['end1'] - result['start1'])
    assert (covered_peak_frac <= 1).all() and (covered_peak_frac >= 0).all()

    is_proper = covered_peak_frac == 1
    print(f'Improper interactions frac: {1 - is_proper.mean()}')
    result = result[is_proper]

    result['name2'] = drop_id_version(result['name2'])

    result = summarize_pairwise(result[['name1', 'name2']], symmetrize=False)

    return result
