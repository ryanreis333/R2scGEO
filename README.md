# R2scGEO

Quickly search, download, and load public single-cell RNA-seq datasets from
[NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/), straight into R.

It wraps the NCBI E-utilities API (search + metadata) and the GEO FTP server
(supplementary files), with on-disk caching, and parses 10x-style sparse
matrices into a `dgCMatrix`, a `Seurat` object, or a `SingleCellExperiment`.

## Install

```r
# install.packages("devtools")
devtools::install_local("R2scGEO")   # or install_github("ryanreis333/R2scGEO")
```

Core dependencies: `httr2`, `jsonlite`, `Matrix`, `tibble`, `rlang`.
`Seurat` / `SingleCellExperiment` are only needed if you request those output
types.

## Setup (recommended)

NCBI asks API callers to identify themselves; an API key raises your rate limit
from 3 to 10 requests/sec.

```r
Sys.setenv(ENTREZ_EMAIL = "you@example.com")
Sys.setenv(ENTREZ_KEY   = "your_ncbi_api_key")  # optional
```

## Usage

```r
library(R2scGEO)

# 1. Search for human glioblastoma single-cell studies
hits <- scgeo_search("glioblastoma", organism = "Homo sapiens", max_results = 10)
hits[, c("accession", "title", "n_samples")]

# 2. Inspect one study's metadata
scgeo_metadata("GSE164897")

# 3. See what supplementary files are available
scgeo_suppl_files("GSE164897")

# 4. Download supplementary files (cached on disk)
files <- scgeo_download("GSE164897")

# 5. Load — format is auto-detected, whatever GEO provided
obj <- scgeo_load(dirname(files[1]), as = "seurat")
```

### Universal loading

`scgeo_load()` points at a file *or* a directory (e.g. the whole download
folder) and figures out the format itself:

```r
scgeo_load("GSE164897/")                      # scan a dir, auto-detect
scgeo_load("counts.csv.gz")                   # dense table
scgeo_load("filtered_feature_bc_matrix.h5")   # 10x HDF5
scgeo_load("adata.h5ad", as = "sce")          # AnnData -> SingleCellExperiment
scgeo_load("GSE..._RAW.tar")                  # extracted + scanned recursively
scgeo_load(dir, sample = "GSM2")              # pick one of several samples
```

If a directory or archive holds several samples, you get a **named list** of
objects; pass `sample=` to select just one.

### Supported formats

| Input | Detected as | Reader / backend |
|---|---|---|
| `matrix.mtx` + `barcodes.tsv` + `features/genes.tsv` (± `.gz`) | `mtx` | `Matrix::readMM` |
| `*.h5` (CellRanger) | `h5` | `Seurat::Read10X_h5`, else `hdf5r` |
| `*.h5ad` (AnnData) | `h5ad` | `zellkonverter` → `anndata` → manual `hdf5r` |
| `*.loom` | `loom` | `hdf5r` |
| `*.rds` (Seurat / SCE / matrix) | `rds` | `readRDS` + counts extraction |
| `*.csv` / `*.tsv` / `*.txt` (± `.gz`) | `dense` | `data.table::fread`, else `read.delim` |
| `*.tar` / `*.tgz` / `*.zip` | archive | extracted, then re-scanned recursively |
| extensionless delimited text | `dense` | content sniffing |

All readers normalize to a genes × cells sparse `dgCMatrix` before converting
to your requested output type. Optional backends (`hdf5r`, `zellkonverter`,
`anndata`, `data.table`, `Seurat`, `SingleCellExperiment`) are only required for
the formats that use them — install them as needed.

Caching: downloads live in `scgeo_cache_dir()` (override with the `SCGEO_CACHE`
env var). Clear with `scgeo_clear_cache()` or `scgeo_clear_cache("GSE164897")`.

## Functions

| Function | Purpose |
|---|---|
| `scgeo_search()` | Search GEO series; returns a metadata tibble |
| `scgeo_metadata()` | Metadata for specific accessions |
| `scgeo_suppl_files()` | List an accession's supplementary files |
| `scgeo_download()` | Download (and cache) supplementary files |
| `scgeo_load()` | Auto-detect any format and load to matrix / Seurat / SCE |
| `scgeo_load_10x()` | Load an mtx triplet (thin wrapper, back-compat) |
| `scgeo_cache_dir()`, `scgeo_clear_cache()` | Manage the cache |

## Development & testing

```r
# from the package directory
devtools::document()   # generate man/*.Rd from roxygen comments (first run)
devtools::test()       # run the testthat suite (offline; uses fixtures + mocks)
devtools::check()      # full R CMD check
```

Or from a shell:

```sh
R -e 'devtools::document()'
R CMD build R2scGEO
R CMD check R2scGEO_0.1.0.tar.gz
```

The test suite is fully offline: network calls are mocked and the loader is
tested against a real on-disk Matrix Market triplet, so `devtools::test()` runs
without hitting NCBI.

## Notes

GEO is heterogeneous — not every series ships a clean 10x triplet. Some provide
a single `*_RAW.tar`, per-sample matrices, or processed `.h5`/`.rds`/`.csv`
files. Use `scgeo_suppl_files()` to inspect before downloading; `scgeo_load_10x()`
targets the common `matrix.mtx` + `barcodes.tsv` + `features.tsv`/`genes.tsv`
layout.
