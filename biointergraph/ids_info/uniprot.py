from re import escape
from urllib.parse import urlencode

import pandas as pd

from ..shared import _read_tsv, memory
from ..ids_mapping import yapid2ids_by_type


GO_NUCLEAR = [
    'nucleus [GO:0005634]',
    'nucleoplasm [GO:0005654]',
    'nucleolus [GO:0005730]',
    'nuclear speck [GO:0016607]',
    'nuclear body [GO:0016604]',
    'nuclear matrix [GO:0016363]',
    'nucleosome [GO:0000786]'
]


def _is_nuclear(data: pd.DataFrame) -> pd.Series:
    data = data[['Gene Ontology (cellular component)', 'Subcellular location [CC]']]
    data = data.fillna('')

    result = data['Subcellular location [CC]'].str.contains('Nucleus') / 2
    regex = '|'.join(map(escape, GO_NUCLEAR))
    result += data['Gene Ontology (cellular component)'].str.contains(regex) / 2

    assert result.notna().all() and result.between(0, 1).all()

    return result


@memory.cache
def uniprot_id_info(organism_id: str = '9606') -> pd.DataFrame:
    params = {
        'query': ' AND '.join(f'{key}:{value}' for key, value in [
            ('reviewed', 'true'),
            ('organism_id', organism_id)
        ]),
        'format': 'tsv',
        'compressed': 'true',
        'fields': ','.join([
            'accession', 'id',
            'ft_dna_bind', 'reviewed',
            'cc_interaction', 'go_c',
            'ft_intramem', 'cc_subcellular_location',
            'ft_topo_dom', 'ft_transmem', 'ft_zn_fing'
        ])
    }

    result = _read_tsv(
        f'https://rest.uniprot.org/uniprotkb/stream?{urlencode(params)}',
        compression='gzip'
    )

    result['is_nuclear'] = _is_nuclear(result)
    return result


def yapid2is_nuclear() -> pd.Series:
    uniprot2is_nuclear = uniprot_id_info()
    uniprot2is_nuclear = uniprot2is_nuclear.set_index('Entry', verify_integrity=True)['is_nuclear']

    result = yapid2ids_by_type()
    result = result['uniprot'].explode()

    result = result.map(uniprot2is_nuclear).fillna(0)
    result = result.groupby(level=0).sum().ge(1)

    return result
