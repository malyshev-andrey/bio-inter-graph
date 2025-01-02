import pandas as pd
import networkx as nx

from ..shared import memory


@memory.cache
def _build_yagid_graph():
    result = pd.Series(['a', 'b', 'c'])
    return result

def id2yagid(ids: pd.Series) -> pd.Series:
    mapping = _build_yagid_graph()
    ids = ids.map(mapping)
    return ids
