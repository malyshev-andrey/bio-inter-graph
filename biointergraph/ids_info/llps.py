import pandas as pd

from ..shared import memory, remote_file2local
from ..ids_mapping import yapid2ids_by_type


AICAP_URL = (
    'https://static-content.springer.com/esm/'
    'art%3A10.1186%2Fs13059-021-02456-2/MediaObjects/'
    '13059_2021_2456_MOESM2_ESM.xlsx'
)
# The workbook also holds 'condition 1(1,6-HD-1)' (fewer proteins) and the
# 'condition 2(2,5-HD)' specificity control; the authors state that every
# figure in the paper is based on condition 2 of the 1,6-HD experiment.
AICAP_SHEET = 'condition 2(1,6-HD-2)'

AICAP_THRESHOLD = 1.0


@memory.cache
def aicap_info() -> pd.DataFrame:
    """
    Load AICAP values for chromatin-associated proteins in K562.

    AICAP (anti-1,6-hexanediol capacity index) is measured by Hi-MS: chromatin-associated
    proteins are extracted before and after 1,6-hexanediol treatment and quantified by mass
    spectrometry. 1,6-hexanediol dissolves the weak multivalent interactions that hold
    biomolecular condensates together, so a protein retained on chromatin through phase
    separation is washed out and gets a LOW AICAP value, while a protein held by a stable
    complex keeps a HIGH one. The scale is therefore inverted with respect to the
    "phase separation propensity" it is used as a proxy for.

    Source: Shi et al., Genome Biology 2021, doi:10.1186/s13059-021-02456-2 (CC BY 4.0),
    Additional file 2. Covers only the chromatin-associated fraction of the proteome.

    Returns
    -------
    pd.DataFrame
        Indexed by UniProt accession, with columns 'aicap' (float) and 'pvalue'
        (t-test p-value reported by the authors).
    """
    path = remote_file2local(AICAP_URL).removeprefix('file://')

    result = pd.read_excel(path, sheet_name=AICAP_SHEET, dtype='str')

    result = result[['Uniprot', 'AICAP', 'AICAP t-test Pvalue']]
    result.columns = ['uniprot', 'aicap', 'pvalue']

    result = result.dropna(subset=['uniprot', 'aicap'])

    # a few accessions in the published table carry stray tabs
    result['uniprot'] = result['uniprot'].str.strip()

    regex = r'^([A-Z0-9]{6}|[A-Z0-9]{10})$'
    assert result['uniprot'].str.match(regex).all()

    result = result.astype({'aicap': 'float', 'pvalue': 'float'})
    assert (result['aicap'] > 0).all()

    result = result.sort_values('aicap')
    result = result.drop_duplicates('uniprot', keep='first')
    result = result.set_index('uniprot', verify_integrity=True)

    return result


def yapid2aicap(ids: pd.Series|None = None) -> pd.Series:
    """
    Map AICAP values onto YAPID nodes.

    A YAPID may merge several UniProt accessions; the minimum is kept, consistently with
    `yapid2is_disordered` taking the maximum disorder fraction: both pick the most
    condensate-like member of the group.

    Proteins outside the assayed chromatin-associated fraction are absent from the result
    rather than filled with a default — a missing AICAP means "not measured", and on this
    inverted scale any fill value would read as an actual phase separation propensity.
    """
    uniprot2aicap = aicap_info()['aicap']

    result = yapid2ids_by_type()['uniprot'].explode()
    result = result.map(uniprot2aicap)
    result = result.groupby(level=0).min()
    result = result.dropna()
    result.name = 'aicap'

    if ids is not None:
        result = ids.map(result)

    return result


def yapid2is_llps(
        ids: pd.Series|None = None, *,
        threshold: float = AICAP_THRESHOLD
    ) -> pd.Series:
    """
    Dichotomize AICAP into a phase separation flag, using the cutoff of the original paper.

    Returns a nullable boolean Series: pd.NA marks proteins with no AICAP measurement,
    which must not be conflated with proteins measured and found insensitive to 1,6-hexanediol.
    """
    result = yapid2aicap()
    result = result.lt(threshold).astype('boolean')
    result.name = 'is_llps'

    if ids is not None:
        result = ids.map(result).astype('boolean')

    return result
