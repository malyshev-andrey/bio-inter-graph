import pandas as pd
import networkx as nx

from .ENCODE import encode_eCLIP2pairwise
from .karr_seq import load_karr_seq_data
from .ric_seq import load_ric_seq_data
from .protein import load_IntAct_interactions, load_biogrid_interactions, load_string_interactions
from ..shared import memory


@memory.cache
def build_main_graph():
    data = [
        encode_eCLIP2pairwise(assembly='hg38', annotation='gencode', pvalue=0.05),
        load_ric_seq_data(0.05),
        load_karr_seq_data(0.05),
        load_IntAct_interactions(),
        load_biogrid_interactions(),
        load_string_interactions()
    ]

    for df in data:
        assert df.shape[1] == 2
        df.columns = 'source', 'target'

    data = pd.concat(data)
    assert data.shape[1] == 2
    regex = r'^YA[PG]ID\d{7}$'
    assert (
        data['source'].str.match(regex).all() and
        data['target'].str.match(regex).all()
    )

    graph = nx.from_pandas_edgelist(data)
    return graph
