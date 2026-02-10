import pandas as pd

from ..shared import memory, _read_tsv
from ..ids_mapping import id2yapid


def _to_pairwise(id1: pd.Series, id2: pd.Series) -> pd.DataFrame:
    result = pd.DataFrame({
        'yapid1': id2yapid(id1),
        'yapid2': id2yapid(id2)
    })
    assert result['yapid1'].str.startswith('YAPID').all()
    assert result['yapid2'].str.startswith('YAPID').all()

    result = pd.concat([
        result,
        result.rename(columns={'yapid1': 'yapid2', 'yapid2': 'yapid1'})
    ])
    result = result[result['yapid1'] < result['yapid2']]
    result = result.drop_duplicates()
    return result


@memory.cache
def load_biogrid_interactions() -> pd.DataFrame:

    def _is_appropriate_qualifications(qual: pd.Series) -> pd.Series:
        qual = qual.str.lower()
        qual = qual.str.replace(r'[^a-z0-9]', ' ', regex=True)
        result = ~qual.str.contains(r'\b(?:proximity ligation|pla|chip)\b', regex=True)
        return result

    result = _read_tsv(
        'https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-MV-Physical-LATEST.tab3.zip',
        filter_func=lambda df: df[
            df['Organism ID Interactor A'].eq('9606') &
            df['Organism ID Interactor B'].eq('9606') &
            ~df['Experimental System'].isin({'Protein-RNA', 'Affinity Capture-RNA'}) &
            _is_appropriate_qualifications(df['Qualifications'])
        ]
    )

    assert result['Organism Name Interactor A'].eq('Homo sapiens').all()
    assert result['Organism Name Interactor B'].eq('Homo sapiens').all()

    assert not result['Experimental System'].str.contains('RNA').any()
    assert result['Experimental System Type'].eq('physical').all()

    result = _to_pairwise(
        result['BioGRID ID Interactor A'],
        result['BioGRID ID Interactor B']
    )

    return result


@memory.cache
def load_string_interactions(min_score: int = 700) -> pd.DataFrame:
    result = _read_tsv(
        f'https://stringdb-downloads.org/download/stream/protein.physical.links.v12.0/9606.protein.physical.links.v12.0.min{min_score}.onlyAB.tsv.gz',
        usecols=['protein1', 'protein2']
    )

    id_regex = r'^9606.ENSP\d{11}$'
    assert result['protein1'].str.match(id_regex).all()
    assert result['protein2'].str.match(id_regex).all()

    result['protein1'] = result['protein1'].str.removeprefix('9606.')
    result['protein2'] = result['protein2'].str.removeprefix('9606.')

    result = _to_pairwise(result['protein1'], result['protein2'])
    return result


@memory.cache
def load_intact_interactions() -> pd.DataFrame:
    result = _read_tsv(
        'https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/species/human.txt',
        usecols=[
            '#ID(s) interactor A',
            'ID(s) interactor B',
            'Interaction detection method(s)',
            'Publication Identifier(s)',
            'Taxid interactor A',
            'Taxid interactor B',
            'Interaction type(s)',
            'Type(s) interactor A',
            'Type(s) interactor B',
            'Confidence value(s)'
        ],
        filter_func=lambda df: df[
            df['Type(s) interactor A'].eq('psi-mi:"MI:0326"(protein)') &
            df['Type(s) interactor B'].eq('psi-mi:"MI:0326"(protein)') &
            df['Taxid interactor A'].eq('taxid:9606(human)|taxid:9606(Homo sapiens)') &
            df['Taxid interactor B'].eq('taxid:9606(human)|taxid:9606(Homo sapiens)') &
            df['Interaction type(s)'].isin({
                'psi-mi:"MI:0914"(association)',
                'psi-mi:"MI:0915"(physical association)',
                'psi-mi:"MI:0407"(direct interaction)'
            })
        ],
        na_values=['-'],
        use_cache=True
    )

    result['PMID'] = result['Publication Identifier(s)'].str.extract(r'pubmed:(\d+)', expand=False)
    result = result[~result['PMID'].isna()]
    assert not result['Interaction detection method(s)'].isna().any()

    for c in '#ID(s) interactor A', 'ID(s) interactor B':
        result = result[result[c].str.startswith('uniprotkb:')]
        result[c] = result[c].str.removeprefix('uniprotkb:')
        result[c] = result[c].str.split('-', expand=True)[0]
        assert result[c].str.match(r'^([A-Z0-9]{6}|[A-Z0-9]{10})$').all()

    result['yapid1'] = id2yapid(result['#ID(s) interactor A'])
    result['yapid2'] = id2yapid(result['ID(s) interactor B'])
    result = result[
        result['yapid1'].str.startswith('YAPID') &
        result['yapid2'].str.startswith('YAPID')
    ]

    result = pd.concat([
        result,
        result.rename(columns={'yapid1': 'yapid2', 'yapid2': 'yapid1'})
    ])
    result = result[result['yapid1'] < result['yapid2']]

    result['weight'] = result['Confidence value(s)'].str.extract(r'intact\-miscore:([0-9\.]+)').astype('float')

    result = result.grouby(
        ['yapid1', 'yapid2', 'PMID', 'Interaction detection method(s)'],
        as_index=False,
        observed=True
    )['weight'].max()
    result = result.groupby(['yapid1', 'yapid2'], as_index=False, observed=True).agg(
        size=('weight', 'size'),
        weight=('weight', 'max')
    )
    result = result.loc[result['size'] > 1, ['yapid1', 'yapid2', 'weight']]

    return result
