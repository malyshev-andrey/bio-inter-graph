# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`bio-inter-graph` is a data-integration library that builds one heterogeneous interaction graph for the human K562 cell line out of ~15 public datasets. Nodes belong to three namespaces, minted by this package rather than taken from any external database:

| Prefix | Entity | Minted from |
|---|---|---|
| `YAGID#######` | gene/RNA | connected components of a cross-reference graph over Entrez, ENSG, ENST, RefSeq NM/NR and extended-annotation `EXTG` ids |
| `YAPID#######` | protein | connected components over UniProt accessions, ENSP, BioGRID ids and `SYMBOL:<gene name>` |
| `YALID#######` | DNA locus | 500 bp bins of the K562 ChromHMM segmentation (18 states collapsed to 6, ENCODE blacklist removed, SPIN state attached) |

Every interaction loader ends by mapping its native identifiers into these namespaces, which is what makes edges from RNA-RNA, RNA-protein, protein-DNA and DNA-DNA assays composable into a single graph.

## Commands

```bash
# environment actually used on this machine (uv-created; system python3 has no networkx)
~/biointergraph/.venv/bin/python -c "import biointergraph; print(biointergraph.__version__)"

# tests (pytest is NOT installed in that venv yet)
pip install -e '.[dev]'
pytest
pytest tests/test_shared_remote_file2local.py::test_remote_file2local_downloads_and_caches
```

No linter, formatter or type checker is configured; there is no CI. The test suite covers only `remote_file2local`.

`__version__` comes from installed distribution metadata (`importlib.metadata`), so the package must be installed for import to work at all. The venv holds a **non-editable** install pinned to an older commit — running with the repo as cwd shadows it with the working tree, running from anywhere else silently uses the stale snapshot.

## Architecture

### Layering (import direction is strict)

```
shared.py  →  ids/  →  annotations/  →  ids_mapping/  →  ids_info/  →  interactions/
```

Modules import only from layers to their left; `ids_info.main` needs `ids_mapping`, `ids_mapping` needs `annotations`, and `interactions` needs everything. Adding a leftward import creates a circular import (one has already been fixed once) — if a low layer needs something from a high one, import it inside the function.

### The loader contract

Every `load_*_data` / `load_*_interactions` function in `interactions/` returns a DataFrame with **exactly three columns**: two id columns (`YA?ID` values) and `weight`, one row per unordered pair. `graph.py::_wrapper` enforces this with assertions, then:

1. renames the id columns to `source`/`target` and sorts each pair lexicographically (`source < target`);
2. writes the per-dataset edge list to `~/.cache/bio-inter-graph/datasets/<dataset>-<args>.tsv.gz`;
3. **replaces `weight` by its quantile rank** `rank/(n+1) ∈ (0,1)` within that dataset.

Raw weights are heterogeneous by construction (`-log10(p)` for KARR-seq / RIC-seq / PRIM-seq / Hi-C, peak scores for CLIP-like assays, a blended `log1p` of experiment counts for GTRD, a constant for Red-C), so the quantile ranking is what makes them comparable. After dedup across datasets the edge keeps `max` of the ranks and a comma-joined `dataset` string.

Consequence for any analysis: `weight` is a **similarity/strength** in (0,1), not a distance. NetworkX shortest-path functions (`betweenness_centrality`, `closeness_centrality`) treat `weight` as edge *length*, so convert first (e.g. `1 - weight`).

### From genomic intervals to node ids

Interval-based assays go through `interactions/main.py::_annotate_peaks` → `annotations/intersect.py::best_left_intersect`, which is a PyRanges left join keeping the single best overlap per left interval, ranked by Jaccard index, with unmatched rows dropped. Chromosome naming is unified through the UCSC `chromAlias` table (`annotations/ucsc.py::unify_chr`). Pair-count significance for cluster/chimera data comes from `interactions/main.py::summarize_pairwise` (Fisher exact + BH correction + PMI).

### Caching (two independent layers)

- `@memory.cache` — joblib `Memory` in `~/.cache/bio-inter-graph/joblib`, keyed on function arguments. Applied to almost every expensive loader, so a signature change silently invalidates a very expensive cache entry.
- `remote_file2local()` / `_read_tsv(use_cache=True)` — fsspec `simplecache` in `~/.cache/bio-inter-graph/fsspec`, keyed on canonicalized URL, writing a `.meta.json` next to each download. Handles `::`-chained fsspec URLs. `_clear_fsspec_cache()` wipes it.

A populated cache is large: on this machine 17 GB of fsspec downloads, 1.3 GB joblib, 315 MB per-dataset edge lists.

### Graph construction

`build_main_graph()` and `build_light_graph()` **do not rebuild by default** — they fetch a prebuilt `edges.tsv.gz` / `edges_light.tsv.gz` asset from the latest GitHub release of `malyshev-andrey/bio-inter-graph`. `rebuild=True` re-downloads and recomputes every dataset (hours, tens of GB) and writes the edge list to `~/.cache/bio-inter-graph/`, from where it is uploaded as a release asset manually.

Both builds assert a **single connected component** (`_remove_minor_components` drops everything but the largest). `build_light_graph` additionally removes mRNA nodes, degree-1 DNA nodes, and nodes with no RNA/protein neighbour.

Scale (measured 2026-09-10): main graph 4 196 699 nodes / 31 932 731 edges; light graph 2 078 151 / 16 001 871. Node counts YALID / YAGID / YAPID = 4 091 070 / 87 788 / 17 841. Edges are dominated by Hi-C DNA-DNA pairs (56 %) and ChIP-seq DNA-protein pairs (24 %); degree is wildly skewed (YAPID mean 507, max 261 838; YALID mean 11.6). Raw degree therefore mostly reflects **which entities were assayed**, not biology — normalise within node type or per neighbour-type channel (`graph.py::_node2neighbors_types`) before calling anything a hub. Exact betweenness is infeasible at this size (O(V·E) ≈ 10¹⁴).

### Precomputed static resources

`biointergraph/static/` holds the artefacts that define node identity: `id2yagid.json`, `id2yapid.json`, `chromhmm_500.tsv.gz`. The functions that would rebuild them are gated behind hardcoded `REBUILD_YAGID_MAPPING` / `REBUILD_YAPID_MAPPING` / `REBUILD_CHROMHMM_ANNOTATION = False` constants, and the ChromHMM table is additionally checked against `CHROMHMM_500_HASH`. Regenerating any of them renumbers node ids, invalidating every cached edge list, the published release assets and any saved analysis output.

### Analysis surface

`interactions/graph.py` and `interactions/analysis.py` provide `describe_graph`, `describe_nodes` (attaches biotype, ChromHMM state, neighbour-type counts), `describe_edges`, `detect_communities` (Louvain, optional g:Profiler GO enrichment and nuclear fraction), `graph2random_walks`, `indirect_interactions` (common-neighbour features for non-adjacent pairs), and `graph_datasets_stats(latex='en'|'ru')` which renders the dataset provenance table for the write-up.

Node attribute providers live in `ids_info/`: `yagid2biotype` (weighted vote across id types), `yapid2is_nuclear` (UniProt GO/subcellular), `yapid2is_disordered` (MobiDB), `yagid2rna_localization` (ENCODE K562 nucleus/cytoplasm RSEM, or APEX-seq). The last two are implemented but not yet consumed anywhere.

Exploratory analysis is not kept in this repo: standalone scripts (`analyze_*.py` with argparse writing `metadata.json` + CSVs + a log) and notebooks live in dated folders under `~/biointergraph/`.

## Conventions and traps

- `_read_tsv` reads everything as `dtype='str'` and chunks with a tqdm progress bar by default; pass `chunksize=None` for whole-file reads, and cast numeric columns explicitly.
- `assert` is used as the data-validation mechanism throughout (id regexes, coordinate sanity, join cardinality, expected state names). An assertion failure normally means an upstream source changed its schema, URL or contents — investigate the source rather than relaxing the assertion.
- Data source URLs are hardcoded, including Google Drive file ids (`GOOGLE_DRIVE_URL`), GEO/FTP paths and `http://gtrd.biouml.org:8888`. `gtrd.py` downloads the UCSC `bigBedToBed` binary at runtime (Linux x86_64 only) and `prim_seq.py` fetches over HTTPS with `verify=False`.
- Human/K562 is hardcoded in places (UniProt/MobiDB `organism_id='9606'`, ChromHMM sample `BSS00762`).
- `load_rdsprite_data` (`interactions/sprite.py`) is exported but is **not** a graph loader: it returns cluster-level rows with tuple columns, not pairwise weighted edges, and is absent from the dataset list in `build_main_graph`. Its source sample is H1, not K562.
- When adding an interaction dataset: new module in `interactions/`, `@memory.cache`d entry point returning the 3-column frame, export it in `interactions/__init__.py`, add it to the `data` list in `build_main_graph`, and add a row to the `metadata` dict in `graph_datasets_stats` (otherwise it appears there with empty provenance).
- `requests-cache` is still declared as a dependency but is no longer imported anywhere; `ids_mapping/uniprot.py` is an empty placeholder.
