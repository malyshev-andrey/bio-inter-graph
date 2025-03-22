import tempfile
import zipfile
from typing import Callable

import requests
import pandas as pd
from tqdm.auto import tqdm

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

    assert (
        result['Organism Name Interactor A'].eq('Homo sapiens') &
        result['Organism Name Interactor B'].eq('Homo sapiens')
    ).all()
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
def load_IntAct_interactions(
        filter_func: Callable[[pd.DataFrame], pd.DataFrame]|None = None
    ) -> pd.DataFrame:
    if filter_func is None:
        filter_func = lambda df: df[
            df['Type(s) interactor A'].eq('psi-mi:"MI:0326"(protein)') &
            df['Type(s) interactor B'].eq('psi-mi:"MI:0326"(protein)') &
            df['ID(s) interactor A'].str.startswith('uniprotkb:') &
            df['ID(s) interactor B'].str.startswith('uniprotkb:') &
            df['Taxid interactor A'].eq('taxid:9606(human)|taxid:9606(Homo sapiens)') &
            df['Taxid interactor B'].eq('taxid:9606(human)|taxid:9606(Homo sapiens)') &
            df['Interaction type(s)'].isin([
                'psi-mi:"MI:0914"(association)',
                'psi-mi:"MI:0915"(physical association)',
                'psi-mi:"MI:0407"(direct interaction)'
            ]) &
            (df['ID(s) interactor A'] != df['ID(s) interactor B'])
        ]

    url = 'https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/species/human.zip'
    with tempfile.NamedTemporaryFile(suffix=".zip", delete=True) as db_archive:
        response = requests.get(url, stream=True)
        response.raise_for_status()

        total_size = int(response.headers.get("content-length", 0))
        with tqdm(desc=url, total=total_size, unit="B", unit_scale=True) as progress_bar:
            for chunk in response.iter_content(chunk_size=8192):
                progress_bar.update(len(chunk))
                db_archive.write(chunk)
        db_archive.flush()

        with zipfile.ZipFile(db_archive.name, 'r') as z:
            with z.open('human.txt') as db:
                result = pd.read_csv(db, sep='\t')
    result = result.rename(columns={'#ID(s) interactor A': 'ID(s) interactor A'})
    result = filter_func(result)

    uniprot_id_regex = r'^uniprotkb:(.{6}(-(\d+|PRO_\d{10}))?|.{10})$'
    assert (
        result['ID(s) interactor A'].str.match(uniprot_id_regex).all() and
        result['ID(s) interactor B'].str.match(uniprot_id_regex).all()
    )
    result['ID(s) interactor A'] = result['ID(s) interactor A'].str.removeprefix('uniprotkb:')
    result['ID(s) interactor B'] = result['ID(s) interactor B'].str.removeprefix('uniprotkb:')

    result = result[['ID(s) interactor A', 'ID(s) interactor B', 'Publication Identifier(s)']]

    n = result.shape[0]
    ids = ['ID(s) interactor A', 'ID(s) interactor B']
    swap = dict(zip(ids, ids[::-1]))
    result = pd.concat([
        result,
        result.rename(columns=swap)
    ])
    result = result[result['ID(s) interactor A'] < result['ID(s) interactor B']]
    assert result.shape[0] == n

    result['PMID'] = result['Publication Identifier(s)'].str.extract('pubmed:(\d+)')
    result = result[~result['PMID'].isna()]
    result['n_publications'] = result.groupby(ids)['PMID'].transform('nunique')

    result = result[result['n_publications'] > 1]

    result['yapid1'] = id2yapid(result['ID(s) interactor A'])
    result['yapid2'] = id2yapid(result['ID(s) interactor B'])
    result = result[
        result['yapid1'].str.startswith('YAPID') &
        result['yapid2'].str.startswith('YAPID')
    ]

    ids = ['yapid1', 'yapid2']
    result = result[ids]
    swap = dict(zip(ids, ids[::-1]))
    result = pd.concat([
        result,
        result.rename(columns=swap)
    ])
    result = result[result['yapid1'] < result['yapid2']]
    result = result.drop_duplicates()
    return result
