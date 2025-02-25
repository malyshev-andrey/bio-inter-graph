import tempfile
import tarfile
import sqlite3

import requests
import pandas as pd

from ..shared import memory


def _load_OrgHsEgDb(query_func):
    url = 'https://bioconductor.org/packages/release/data/annotation/src/contrib/org.Hs.eg.db_3.20.0.tar.gz'
    with tempfile.NamedTemporaryFile(suffix=".tar.gz", delete=True) as db_archive:
        response = requests.get(url, stream=True)
        response.raise_for_status()
        for chunk in response.iter_content(chunk_size=8192):
            db_archive.write(chunk)
        db_archive.flush()

        with tempfile.NamedTemporaryFile(suffix='.sqlite', delete=True) as db_file:
            with tarfile.open(db_archive.name, 'r:gz') as tar:
                db_path = 'org.Hs.eg.db/inst/extdata/org.Hs.eg.sqlite'
                db_file.write(tar.extractfile(db_path).read())

            with sqlite3.connect(db_file.name) as conn:
                return query_func(conn)


@memory.cache
def load_OrgHsEgDb_pairwise(id1_type: str, id2_type: str) -> pd.DataFrame:
    def query_func(conn) -> pd.DataFrame:
        def load_table(id_type: str) -> pd.DataFrame:
            aliases = {
                'gene_id': 'entrezgene_id',
                'accession': 'refseq_transcript_id',
                'ensembl_id': 'ensembl_gene_id',
                'trans_id': 'ensembl_transcript_id'
            }

            result = None
            if id_type == 'entrezgene_id':
                result = pd.read_sql_query("SELECT * FROM genes", conn)
            elif id_type == 'refseq_transcript_id':
                result = pd.read_sql_query("SELECT * FROM refseq", conn)
                result = result[result['accession'].str[:2].isin({'NR', 'NM'})]
            elif id_type == 'ensembl_gene_id':
                result = pd.concat([
                    pd.read_sql_query("SELECT * FROM ensembl", conn),
                    pd.read_sql_query("SELECT * FROM ensembl2ncbi", conn),
                    pd.read_sql_query("SELECT * FROM ncbi2ensembl", conn)
                ])
            elif id_type == 'ensembl_transcript_id':
                result = pd.read_sql_query("SELECT * FROM ensembl_trans", conn)
            else:
                assert id_type in aliases.values()

            result = result.drop_duplicates()
            result = result.rename(columns=aliases)

            return result


        result = load_table(id1_type).merge(load_table(id2_type), how='inner')
        result = result.drop(columns=['_id'])
        result = result.drop_duplicates()

        assert result.shape[1] == 2

        return result

    result = _load_OrgHsEgDb(query_func)
    print(f'OrgHsEgDb query ({id1_type}, {id2_type}) result: {result.shape}')
    return result


@memory.cache
def entrezgene_id2go() -> pd.DataFrame:
    def query_func(conn):
        query = """
            SELECT gene_id, go_id, "bp" as category FROM genes INNER JOIN go_bp USING(_id)
        UNION ALL
            SELECT gene_id, go_id, "mf" as category FROM genes INNER JOIN go_mf USING(_id)
        UNION ALL
            SELECT gene_id, go_id, "cc" as category FROM genes INNER JOIN go_cc USING(_id)
        """
        result = pd.read_sql_query(query, conn, dtype='str')

        result = result.drop_duplicates()

        return result

    return _load_OrgHsEgDb(query_func)
