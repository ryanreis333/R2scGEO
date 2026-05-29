# ---------------------------------------------------------------------------
# Format-specific readers. Each returns a genes x cells sparse dgCMatrix
# (or, for .rds, possibly a Seurat/SCE object passed straight through).
# ---------------------------------------------------------------------------

# Coerce anything matrix-like to a named dgCMatrix.
to_dgc <- function(x) {
  if (methods::is(x, "dgCMatrix")) return(x)
  if (is.data.frame(x)) x <- as.matrix(x)
  if (methods::is(x, "Matrix")) {
    x <- methods::as(x, "generalMatrix")
    if (!methods::is(x, "CsparseMatrix")) x <- methods::as(x, "CsparseMatrix")
  } else {
    x <- Matrix::Matrix(x, sparse = TRUE)
    if (!methods::is(x, "CsparseMatrix")) {
      x <- methods::as(methods::as(x, "generalMatrix"), "CsparseMatrix")
    }
  }
  methods::as(x, "dgCMatrix")
}

# --- 10x Matrix Market triplet ---------------------------------------------
read_triplet <- function(trip) {
  if (is.na(trip$matrix)) {
    rlang::abort("No matrix.mtx file found.", class = "scgeo_no_matrix")
  }
  mat <- to_dgc(Matrix::readMM(open_maybe_gz(trip$matrix)))
  feats <- read_lines_col(trip$features, 2, fallback_col = 1)
  barc  <- read_lines_col(trip$barcodes, 1)
  if (length(feats) != nrow(mat)) feats <- paste0("gene", seq_len(nrow(mat)))
  if (length(barc)  != ncol(mat)) barc  <- paste0("cell", seq_len(ncol(mat)))
  rownames(mat) <- make.unique(feats)
  colnames(mat) <- make.unique(barc)
  mat
}

# --- Dense delimited table (genes x cells) ---------------------------------
read_dense <- function(path) {
  sep <- guess_sep(path)
  first <- read_table_row(path, sep, skip = 0)
  second <- read_table_row(path, sep, skip = 1)

  if (ncol(second) == ncol(first) + 1) {
    cells <- as.character(first[1, ])
    df <- utils::read.table(open_maybe_gz(path), sep = sep, header = FALSE,
                            skip = 1, quote = "\"'", comment.char = "",
                            check.names = FALSE, stringsAsFactors = FALSE)
    genes <- as.character(df[[1]])
    df <- df[, -1, drop = FALSE]
    if (length(cells) == ncol(df)) names(df) <- make.unique(cells)
  } else {
    df <- utils::read.table(open_maybe_gz(path), sep = sep, header = TRUE,
                            quote = "\"'", comment.char = "",
                            check.names = FALSE, stringsAsFactors = FALSE)
    if (ncol(df) < 2) {
      rlang::abort("Dense table has no data columns.", class = "scgeo_bad_table")
    }
    genes <- as.character(df[[1]])
    df <- df[, -1, drop = FALSE]
  }
  if (ncol(df) < 1) rlang::abort("Dense table has no data columns.", class = "scgeo_bad_table")
  m <- as.matrix(df)
  storage.mode(m) <- "double"
  rownames(m) <- make.unique(genes)
  mat <- to_dgc(m)
  # Heuristic: if it looks transposed (far more rows than typical gene counts
  # and column names look like genes), leave as-is; GEO dense tables are
  # overwhelmingly genes-in-rows, so no auto-transpose by default.
  mat
}

read_table_row <- function(path, sep, skip = 0) {
  utils::read.table(open_maybe_gz(path), sep = sep, header = FALSE,
                    nrows = 1, skip = skip, quote = "\"'",
                    comment.char = "", check.names = FALSE,
                    stringsAsFactors = FALSE)
}

guess_sep <- function(path) {
  con <- open_maybe_gz(path)
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  line <- tryCatch(readLines(con, n = 1, warn = FALSE), error = function(e) "")
  if (grepl("\t", line)) "\t" else if (grepl(",", line)) "," else if (grepl(";", line)) ";" else ""
}

# --- 10x HDF5 (.h5) ---------------------------------------------------------
read_h5 <- function(path) {
  if (requireNamespace("Seurat", quietly = TRUE)) {
    return(to_dgc(Seurat::Read10X_h5(path, use.names = TRUE)))
  }
  need_pkg("hdf5r", context = ".h5 files")
  h5_read_10x_manual(path)
}

# Minimal CellRanger .h5 reader (groups: matrix/{data,indices,indptr,shape,
# barcodes, features/name|genes}).
h5_read_10x_manual <- function(path) {
  f <- hdf5r::H5File$new(path, mode = "r")
  on.exit(f$close_all(), add = TRUE)
  grp <- if (f$exists("matrix")) "matrix" else names(f)[1]
  g <- f[[grp]]
  shape <- as.integer(g[["shape"]]$read())   # (genes, cells), CSC by feature
  mat <- Matrix::sparseMatrix(
    i = as.integer(g[["indices"]]$read()) + 1L,
    p = as.numeric(g[["indptr"]]$read()),
    x = as.numeric(g[["data"]]$read()),
    dims = shape, repr = "C"
  )
  feats <- if (g$exists("features")) {
    fg <- g[["features"]]
    if (fg$exists("name")) fg[["name"]]$read() else fg[["id"]]$read()
  } else if (g$exists("genes")) g[["genes"]]$read() else paste0("gene", seq_len(shape[1]))
  barc <- g[["barcodes"]]$read()
  rownames(mat) <- make.unique(as.character(feats))
  colnames(mat) <- make.unique(as.character(barc))
  to_dgc(mat)
}

# --- AnnData (.h5ad) --------------------------------------------------------
read_h5ad <- function(path) {
  if (requireNamespace("zellkonverter", quietly = TRUE)) {
    sce <- zellkonverter::readH5AD(path)
    return(extract_counts(sce))
  }
  if (requireNamespace("anndata", quietly = TRUE)) {
    ad <- anndata::read_h5ad(path)
    return(to_dgc(Matrix::t(to_dgc(ad$X))))   # AnnData is cells x genes
  }
  need_pkg("hdf5r", context = ".h5ad files")
  h5ad_read_manual(path)
}

# Manual AnnData reader: X (CSR/CSC/dense), var/_index (genes), obs/_index (cells).
h5ad_read_manual <- function(path) {
  f <- hdf5r::H5File$new(path, mode = "r")
  on.exit(f$close_all(), add = TRUE)
  X <- f[["X"]]
  # Build `gc`: genes x cells, regardless of how X is stored.
  if (methods::is(X, "H5Group")) {              # sparse
    enc <- tryCatch(hdf5r::h5attr(X, "encoding-type"), error = function(e) "csr_matrix")
    shape <- as.integer(X[["shape"]]$read())     # (n_obs, n_vars) = (cells, genes)
    n_obs <- shape[1]; n_var <- shape[2]
    data <- as.numeric(X[["data"]]$read())
    indices <- as.integer(X[["indices"]]$read()) + 1L
    indptr <- as.numeric(X[["indptr"]]$read())
    if (grepl("csr", enc)) {
      # CSR over cells: indptr per cell, indices = gene. This is exactly the
      # CSC representation of the genes x cells matrix.
      gc <- Matrix::sparseMatrix(i = indices, p = indptr, x = data,
                                 dims = c(n_var, n_obs), repr = "C")
    } else {
      # CSC over genes: indptr per gene, indices = cell -> cells x genes, then t.
      cg <- Matrix::sparseMatrix(i = indices, p = indptr, x = data,
                                 dims = c(n_obs, n_var), repr = "C")
      gc <- Matrix::t(cg)
    }
  } else {
    gc <- Matrix::t(methods::as(X$read(), "CsparseMatrix"))  # dense cells x genes -> t
  }
  m <- gc
  genes <- h5ad_index(f, "var")
  cells <- h5ad_index(f, "obs")
  if (length(genes) == nrow(m)) rownames(m) <- make.unique(genes)
  if (length(cells) == ncol(m)) colnames(m) <- make.unique(cells)
  to_dgc(m)
}

h5ad_index <- function(f, grp) {
  g <- f[[grp]]
  idx <- tryCatch(hdf5r::h5attr(g, "_index"), error = function(e) "_index")
  if (g$exists(idx)) as.character(g[[idx]]$read())
  else if (g$exists("_index")) as.character(g[["_index"]]$read())
  else character(0)
}

# --- Loom (.loom) -----------------------------------------------------------
read_loom <- function(path) {
  need_pkg("hdf5r", context = ".loom files")
  f <- hdf5r::H5File$new(path, mode = "r")
  on.exit(f$close_all(), add = TRUE)
  m <- to_dgc(f[["matrix"]]$read())               # loom matrix is genes x cells
  genes <- loom_attr(f, "row_attrs", c("Gene", "var_names", "gene_names", "Accession"))
  cells <- loom_attr(f, "col_attrs", c("CellID", "obs_names", "cell_names", "Barcode"))
  if (length(genes) == nrow(m)) rownames(m) <- make.unique(genes)
  if (length(cells) == ncol(m)) colnames(m) <- make.unique(cells)
  m
}

loom_attr <- function(f, grp, candidates) {
  if (!f$exists(grp)) return(character(0))
  g <- f[[grp]]
  for (nm in candidates) if (g$exists(nm)) return(as.character(g[[nm]]$read()))
  character(0)
}

# --- Serialized R object (.rds) --------------------------------------------
read_rds <- function(path) readRDS(path)

# Pull a counts matrix out of common single-cell container objects.
extract_counts <- function(obj) {
  if (methods::is(obj, "Matrix") || is.matrix(obj) || is.data.frame(obj)) {
    return(to_dgc(obj))
  }
  if (methods::is(obj, "Seurat")) {
    need_pkg("Seurat")
    return(to_dgc(Seurat::GetAssayData(obj, layer = "counts")))
  }
  if (methods::is(obj, "SingleCellExperiment") || methods::is(obj, "SummarizedExperiment")) {
    need_pkg("SummarizedExperiment")
    a <- SummarizedExperiment::assayNames(obj)
    which <- if ("counts" %in% a) "counts" else a[1]
    return(to_dgc(SummarizedExperiment::assay(obj, which)))
  }
  rlang::abort(sprintf("Don't know how to extract counts from a '%s'.",
                       class(obj)[1]), class = "scgeo_unknown_object")
}
