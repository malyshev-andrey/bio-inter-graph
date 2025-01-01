from urllib.parse import urlencode

import pandas as pd


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
    metadata = pd.read_csv(url, sep='\t', skiprows=1, dtype='str')

    metadata = metadata.loc[:, ~metadata.isna().all()]

    assert metadata.shape[0] > 0, 'No files found!'

    return metadata


def load_encode_eCLIP(**kwargs) -> pd.DataFrame:
    default_kwargs = dict(
        assay='eCLIP',
        processed='true',
        file_format='bed'
    )
    default_kwargs.update(kwargs)
    metadata = _load_encode_metadata(**default_kwargs)

    assert metadata['Biological replicates'].isin({'1', '2', '1,2'}).all()
    assert metadata['Biological replicates'].value_counts(normalize=True).eq(1/3).all()
    metadata = metadata[metadata['Biological replicates'].eq('1,2')]

    result = {}
    for _, row in tqdm(metadata.iterrows(), total=metadata.shape[0]):
        assembly = row['Genome assembly']
        if assembly not in result:
            result[assembly] = []

        url = f'https://www.encodeproject.org{row["Download URL"]}'
        bed = pd.read_csv(
            url, sep='\t', usecols=range(6), header=None,
            names=['chr', 'start', 'end', 'name', 'score', 'strand']
        )
        bed['name'] = row['Target label']
        bed['cell_line'] = row['Biosample name']
        result[assembly].append(bed)

    for assembly in result:
        result[assembly] = pd.concat(result[assembly])
    return result
