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
    return id.str.split('.', expand=True)[0]
