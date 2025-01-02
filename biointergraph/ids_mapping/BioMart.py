from requests.utils import quote
import pandas as pd

from ..shared import ID_TYPES, memory


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

    print(f'BioMart query ({id1_type}, {id2_type}): https://{host}/biomart/martservice?query=...')

    result = pd.read_csv(
        url, sep='\t',
        names=[id1_type, id2_type],
        header=None, dtype='str'
    )

    result = result[~result.isna().any(axis=1)]

    assert not result.duplicated().any()

    return result


@memory.cache
def load_BioMart_pairwise(
        id1_type: str, id2_type: str, *,
        assembly: str = 'GRCh38'
    ) -> pd.DataFrame:

    if id1_type not in ID_TYPES:
        raise ValueError(
            f'"{id1_type}" is not a valid argument. '
            f'Valid arguments are: {", ".join(ID_TYPES)}'
        )
    if id2_type not in ID_TYPES:
        raise ValueError(
            f'"{id2_type}" is not a valid argument. '
            f'Valid arguments are: {", ".join(ID_TYPES)}'
        )
    if id1_type == id2_type:
        raise ValueError(
            f'Expected different id types, got {id1_type} == {id2_type}'
        )

    refseq_map = {
        'refseq_ncrna': 'refseq_transcript_id',
        'refseq_mrna': 'refseq_transcript_id'
    }

    if id1_type == 'refseq_transcript_id':
        refseq_ncrna = _load_BioMart_pairwise('refseq_ncrna', id2_type, assembly=assembly)
        refseq_ncrna = refseq_ncrna.rename(columns=refseq_map)

        refseq_mrna = _load_BioMart_pairwise('refseq_mrna', id2_type, assembly=assembly)
        refseq_mrna = refseq_mrna.rename(columns=refseq_map)

        result = pd.concat([refseq_ncrna, refseq_mrna]).drop_duplicates()
        is_valid = result[id1_type].str[:2].isin({'NM', 'NR'})
        print(f'BioMart: invalid RefSeq transcript ids: {(~is_valid).sum()}')
        result = result[is_valid]
        return result

    elif id2_type == 'refseq_transcript_id':
        refseq_ncrna = _load_BioMart_pairwise(id1_type, 'refseq_ncrna', assembly=assembly)
        refseq_ncrna = refseq_ncrna.rename(columns=refseq_map)

        refseq_mrna = _load_BioMart_pairwise(id1_type, 'refseq_mrna', assembly=assembly)
        refseq_mrna = refseq_mrna.rename(columns=refseq_map)

        result = pd.concat([refseq_ncrna, refseq_mrna]).drop_duplicates()
        is_valid = result[id2_type].str[:2].isin({'NM', 'NR'})
        print(f'BioMart: invalid RefSeq transcript ids: {(~is_valid).sum()}')
        result = result[is_valid]
        return result

    else:
        return _load_BioMart_pairwise(id1_type, id2_type, assembly=assembly)
