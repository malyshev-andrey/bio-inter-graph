from urllib.parse import urlencode

import pandas as pd

from ..shared import _read_tsv, memory


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
    return result
