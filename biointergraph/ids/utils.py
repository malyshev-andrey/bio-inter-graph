import pandas as pd


def is_valid_id(id: pd.Series, id_type: str) -> pd.Series:

    id_type2regexp = {
        'entrezgene_id': r'^\d+$',
        'ensembl_gene_id': r'^ENSG\d{11}$',
        'ensembl_transcript_id': r'^ENST\d{11,12}$',
        'refseq_transcript_id': r'^N[MR]_\d+$'
    }

    if id_type not in id_type2regexp:
        raise ValueError(
            'Expected one of id types:\n' +
            '\n'.join(f'- "{idt}"' for idt in id_type2regexp) + ',\n' +
            f'but got: {id_type}.'
        )

    return id.str.match(id_type2regexp[id_type])


def drop_id_version(id: pd.Series) -> pd.Series:
    """
    Remove version numbers from identifiers in a pandas Series.

    This function processes a pandas Series of identifiers and removes the version number
    (suffix after the first dot `.`) from each identifier. It is typically used for biological
    identifiers like Ensembl Transcript IDs (e.g., ENST00000832827.1) or RefSeq IDs (e.g., NR_024540.1),
    returning the base identifier without the version.

    Parameters
    ----------
    id : pd.Series
        A pandas Series containing identifier strings, some with a version suffix.

    Returns
    -------
    pd.Series
        A pandas Series containing identifiers without version numbers.

    Examples
    --------
    >>> import pandas as pd
    >>> ids = pd.Series(['ENST00000832827.1', 'NR_024540.1', 'ENSG00000198763.2'])
    >>> drop_id_version(ids)
    0    ENST00000832827
    1          NR_024540
    2    ENSG00000198763
    dtype: object
    """
    return id.str.split('.', expand=True)[0]
