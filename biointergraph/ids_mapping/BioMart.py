from requests.utils import quote
import pandas as pd


def _load_BioMart_pairwise(
        id1_type: str, id2_type: str, *,
        assembly: str
    ) -> pd.DataFrame:

    query = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<!DOCTYPE Query>'
        '<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" datasetConfigVersion = "0.6" >'
            '<Dataset name = "hsapiens_gene_ensembl" interface = "default" >'
                f'<Attribute name = "{id1_type}" />'
                f'<Attribute name = "{id2_type}" />'
            '</Dataset>'
        '</Query>'
    )
    host = 'grch37.ensembl.org' if assembly == 'GRCh37' else 'www.ensembl.org'
    url = f'https://{host}/biomart/martservice?query={quote(query)}'

    result = pd.read_csv(
        url, sep='\t',
        names=[id1_type, id2_type],
        header=None, dtype='str'
    )

    result = result[~result.isna().any(axis=1)]

    assert not result.duplicated().any()

    return result


def load_BioMart_pairwise(
        id1_type: str, id2_type: str, *,
        assembly: str = 'GRCh38'
    ) -> pd.DataFrame:
    assert id1_type != id2_type

    refseq_map = {
        'refseq_ncrna': 'refseq_transcript_id',
        'refseq_mrna': 'refseq_transcript_id'
    }

    if id1_type == 'refseq_transcript_id':
        refseq_ncrna = _load_BioMart_pairwise('refseq_ncrna', id2_type, assembly=assembly)
        refseq_ncrna = refseq_ncrna.rename(columns=refseq_map)

        refseq_mrna = _load_BioMart_pairwise('refseq_mrna', id2_type, assembly=assembly)
        refseq_mrna = refseq_mrna.rename(columns=refseq_map)

        return pd.concat([refseq_ncrna, refseq_mrna]).drop_duplicates()

    elif id2_type == 'refseq_transcript_id':
        refseq_ncrna = _load_BioMart_pairwise(id1_type, 'refseq_ncrna', assembly=assembly)
        refseq_ncrna = refseq_ncrna.rename(columns=refseq_map)

        refseq_mrna = _load_BioMart_pairwise(id1_type, 'refseq_mrna', assembly=assembly)
        refseq_mrna = refseq_mrna.rename(columns=refseq_map)

        return pd.concat([refseq_ncrna, refseq_mrna]).drop_duplicates()

    else:
        return _load_BioMart_pairwise(id1_type, id2_type, assembly=assembly)
