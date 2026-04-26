import numpy as np
import pandas as pd

from ..shared import _read_tsv, memory
from ..ids import drop_id_version
from ..ids_mapping import id2yagid


ENCODE_K562_FILES = {
    'nucleus': [
        'https://encode-public.s3.amazonaws.com/2021/04/05/45780b32-7872-47b9-bf87-97436aea87f7/ENCFF501IXI.tsv',
        'https://encode-public.s3.amazonaws.com/2021/04/05/d065de95-8d24-42e2-8048-849977e61f5f/ENCFF142LZQ.tsv',
    ],
    'cytoplasm': [
        'https://encode-public.s3.amazonaws.com/2021/04/08/19b589fe-6676-410d-9a87-39fd1ab7abe4/ENCFF300NAD.tsv',
        'https://encode-public.s3.amazonaws.com/2021/04/08/204f4cc9-6676-4191-ba5e-bc5f039a2821/ENCFF180ZJI.tsv',
    ],
}

APEX_SEQ_DESEQ2_URL = (
    'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116008/suppl/'
    'GSE116008_APEX_seq_deseq2_all_data_polyA_hek293t.txt.gz'
)

RCI_THRESHOLD = 1.0
PSEUDOCOUNT = 0.1
APEX_FDR_THRESHOLD = 0.05


def _load_encode_rsem(url: str) -> pd.Series:
    data = _read_tsv(url, use_cache=True, chunksize=None)
    is_ensg = data['gene_id'].str.startswith('ENSG')
    data.loc[is_ensg, 'gene_id'] = drop_id_version(data.loc[is_ensg, 'gene_id'])
    data['FPKM'] = data['FPKM'].astype(float)
    result = data.groupby('gene_id')['FPKM'].max()
    return result


@memory.cache
def encode_rna_localization() -> pd.Series:
    fpkm = {}
    for fraction, urls in ENCODE_K562_FILES.items():
        replicates = [_load_encode_rsem(url) for url in urls]
        fpkm[fraction] = (replicates[0] + replicates[1]) / 2

    fpkm = pd.DataFrame(fpkm)
    fpkm = fpkm.dropna()

    rci = np.log2(
        (fpkm['cytoplasm'] + PSEUDOCOUNT) /
        (fpkm['nucleus'] + PSEUDOCOUNT)
    )

    result = pd.Series('both', index=rci.index, dtype='object')
    result[rci > RCI_THRESHOLD] = 'cytoplasmic'
    result[rci < -RCI_THRESHOLD] = 'nuclear'

    result = pd.Categorical(result, categories=['nuclear', 'cytoplasmic', 'both'])
    result = pd.Series(result, index=rci.index)
    result.index.name = 'gene_id'

    return result


@memory.cache
def apex_seq_rna_localization() -> pd.Series:
    data = _read_tsv(
        APEX_SEQ_DESEQ2_URL,
        use_cache=True,
        chunksize=None,
        compression='gzip',
    )

    data = data[data['Ensembl_Gene'].str.startswith('ENSG')]
    data = data.set_index('Ensembl_Gene')

    nls_fc = data['NLS.18C'].astype(float)
    nls_p = data['NLS.18C.P'].astype(float)
    nes_fc = data['NES.18C'].astype(float)
    nes_p = data['NES.18C.P'].astype(float)

    nls_sig = (nls_fc > 0) & (nls_p < APEX_FDR_THRESHOLD)
    nes_sig = (nes_fc > 0) & (nes_p < APEX_FDR_THRESHOLD)

    result = pd.Series('both', index=data.index, dtype='object')
    result[nls_sig & ~nes_sig] = 'nuclear'
    result[nes_sig & ~nls_sig] = 'cytoplasmic'

    result = pd.Categorical(result, categories=['nuclear', 'cytoplasmic', 'both'])
    result = pd.Series(result, index=data.index)
    result.index.name = 'gene_id'

    return result


def yagid2rna_localization(source: str = 'encode') -> pd.Series:
    if source == 'encode':
        localization = encode_rna_localization()
    elif source == 'apex':
        localization = apex_seq_rna_localization()
    else:
        raise ValueError(f"source must be 'encode' or 'apex', got '{source}'")

    mapping = id2yagid()

    loc_dict = localization.to_dict()
    ids_with_loc = mapping.index.isin(loc_dict)

    result = pd.DataFrame({
        'yagid': mapping[ids_with_loc].values,
        'loc': [loc_dict[gid] for gid in mapping.index[ids_with_loc]],
    })
    result = result.dropna(subset=['loc'])

    def _aggregate(locs):
        locs = set(locs)
        if len(locs) == 1:
            return locs.pop()
        if 'nuclear' in locs and 'cytoplasmic' in locs:
            return 'both'
        locs.discard('both')
        if len(locs) == 1:
            return locs.pop()
        return 'both'

    result = result.groupby('yagid')['loc'].agg(_aggregate)
    result = result.astype(pd.CategoricalDtype(categories=['nuclear', 'cytoplasmic', 'both']))
    result.name = 'localization'

    return result
