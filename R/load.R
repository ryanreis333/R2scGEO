# ---------------------------------------------------------------------------
# Universal loader: auto-detect and read any GEO scRNA-seq file or directory.
# ---------------------------------------------------------------------------

#' Load scRNA-seq data from any GEO file or directory
#'
#' Auto-detects the storage format and reads it into the object you ask for.
#' Handles 10x Matrix Market triplets, CellRanger `.h5`, AnnData `.h5ad`,
#' `.loom`, serialized `.rds` (Seurat / SingleCellExperiment / matrix), dense
#' delimited tables (`.csv` / `.tsv` / `.txt`, optionally gzipped), and
#' `.tar` / `.zip` archives (which are extracted and scanned recursively).
#'
#' @param path A file or directory (e.g. the output of [scgeo_download()]).
#' @param as Output type: "matrix" (sparse `dgCMatrix`, default), "seurat",
#'   or "sce".
#' @param sample Optional name/regex to pick one triplet when a directory or
#'   archive contains several samples. If omitted and multiple are found, a
#'   named list is returned.
#' @param project Project name passed to Seurat (default "R2scGEO").
#' @return A `dgCMatrix`, `Seurat`, or `SingleCellExperiment` — or a named list
#'   of these when multiple samples are present.
#' @export
#' @examples
#' \dontrun{
#' files <- scgeo_download("GSE164897")
#' obj <- scgeo_load(dirname(files[1]), as = "seurat")
#' }
scgeo_load <- function(path, as = c("matrix", "seurat", "sce"),
                       sample = NULL, project = "R2scGEO") {
  as <- match.arg(as)
  if (!file.exists(path)) rlang::abort(sprintf("Path '%s' not found.", path))

  result <- if (dir.exists(path)) load_dir(path, sample) else load_file(path, sample)

  if (is.list(result) && !methods::is(result, "Matrix")) {
    if (!is.null(sample)) {
      hit <- result[grepl(sample, names(result))]
      if (length(hit) == 0) rlang::abort(sprintf("No sample matched '%s'.", sample))
      result <- hit[[1]]
    } else if (length(result) == 1) {
      result <- result[[1]]
    } else {
      return(lapply(result, as_output, as = as, project = project))
    }
  }
  as_output(result, as, project)
}

#' Load a 10x Matrix Market triplet (kept for back-compatibility)
#'
#' Thin wrapper around [scgeo_load()] for the classic mtx triplet case.
#' @inheritParams scgeo_load
#' @param dir Directory containing the triplet.
#' @return A matrix, Seurat, or SingleCellExperiment object.
#' @export
scgeo_load_10x <- function(dir, as = c("matrix", "seurat", "sce"), project = "R2scGEO") {
  as <- match.arg(as)
  trip <- find_triplets(dir)
  if (length(trip) == 0) {
    rlang::abort("No matrix.mtx file found in directory.", class = "scgeo_no_matrix")
  }
  as_output(read_triplet(trip[[1]]), as, project)
}

# --- dispatch ---------------------------------------------------------------

load_file <- function(path, sample = NULL) {
  fmt <- detect_format(path)
  switch(fmt,
    mtx   = {
      trip <- find_triplets(dirname(path))
      this <- Filter(function(t) normalizePath(t$matrix) == normalizePath(path), trip)
      read_triplet(if (length(this)) this[[1]] else list(matrix = path,
                    barcodes = NA_character_, features = NA_character_))
    },
    h5    = read_h5(path),
    h5ad  = read_h5ad(path),
    loom  = read_loom(path),
    rds   = read_rds(path),
    dense = read_dense(path),
    tar   = ,
    zip   = load_dir(extract_archive(path, fmt), sample),
    rlang::abort(sprintf("Unrecognized file format for '%s'.", basename(path)),
                 class = "scgeo_unknown_format")
  )
}

# Read a directory: prefer triplets; otherwise pick the best single file(s).
load_dir <- function(dir, sample) {
  trips <- find_triplets(dir)
  if (length(trips) > 0) {
    trips <- filter_named_candidates(trips, sample, vapply(trips, function(t) basename(t$matrix), character(1)))
    if (length(trips) == 0) {
      rlang::abort(sprintf("No sample matched '%s'.", sample))
    }
    mats <- lapply(trips, read_triplet)
    names(mats) <- vapply(trips, function(t)
      sub("_?matrix\\.mtx(\\.gz)?$", "", basename(t$matrix), ignore.case = TRUE),
      character(1))
    return(mats)
  }

  fs <- list.files(dir, full.names = TRUE, recursive = TRUE)
  # Recurse into any archives found inside the directory.
  arch <- fs[vapply(fs, function(p) detect_format(p) %in% c("tar", "zip"), logical(1))]
  if (length(arch) > 0) {
    out <- list()
    for (a in arch) out <- c(out, as_list(load_dir(extract_archive(a), sample)))
    if (length(out)) return(out)
  }

  ranked <- rank_data_files(fs)
  ranked <- filter_named_candidates(ranked, sample, basename(ranked))
  if (length(ranked) == 0) {
    msg <- if (is.null(sample)) {
      sprintf("No recognizable data files in '%s'.", dir)
    } else {
      sprintf("No sample matched '%s'.", sample)
    }
    rlang::abort(msg, class = "scgeo_no_files")
  }
  mats <- lapply(ranked, load_file, sample = NULL)
  names(mats) <- basename(ranked)
  mats
}

# Order candidate files by how likely they are to be the count matrix.
rank_data_files <- function(fs) {
  fmt <- vapply(fs, detect_format, character(1))
  keep <- fmt %in% c("mtx", "h5", "h5ad", "loom", "rds", "dense")
  fs <- fs[keep]; fmt <- fmt[keep]
  pri <- c(h5ad = 1, h5 = 2, loom = 3, mtx = 4, rds = 5, dense = 6)
  # Prefer files whose name hints at counts/expression.
  bonus <- ifelse(grepl("count|umi|raw|expr|matrix", basename(fs), ignore.case = TRUE), -0.5, 0)
  fs[order(pri[fmt] + bonus)]
}

as_list <- function(x) if (is.list(x) && !methods::is(x, "Matrix")) x else list(x)

filter_named_candidates <- function(x, sample, labels) {
  if (is.null(sample)) return(x)
  x[grepl(sample, labels)]
}

# --- output coercion --------------------------------------------------------

as_output <- function(x, as, project = "R2scGEO") {
  if (as != "matrix" && (methods::is(x, "Seurat") ||
        methods::is(x, "SingleCellExperiment"))) {
    x <- as_output_passthrough(x, as)
    if (!is.null(x)) return(x)
  }
  mat <- if (methods::is(x, "dgCMatrix")) x else extract_counts(x)
  switch(as,
    matrix = mat,
    seurat = { need_pkg("Seurat"); Seurat::CreateSeuratObject(counts = mat, project = project) },
    sce    = { need_pkg("SingleCellExperiment")
               SingleCellExperiment::SingleCellExperiment(assays = list(counts = mat)) }
  )
}

# If an .rds already holds the requested container type, return it unchanged.
as_output_passthrough <- function(x, as) {
  if (as == "seurat" && methods::is(x, "Seurat")) return(x)
  if (as == "sce" && methods::is(x, "SingleCellExperiment")) return(x)
  NULL
}

# --- shared low-level helpers (used by readers) -----------------------------

open_maybe_gz <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path) else path
}

read_lines_col <- function(path, col, fallback_col = NULL) {
  if (is.na(path) || !file.exists(path)) return(character(0))
  df <- utils::read.delim(open_maybe_gz(path), header = FALSE,
                          stringsAsFactors = FALSE, colClasses = "character")
  if (col > ncol(df) && !is.null(fallback_col)) col <- fallback_col
  if (col > ncol(df)) col <- 1
  df[[col]]
}

need_pkg <- function(pkg, context = "this output type") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    rlang::abort(sprintf("Package '%s' is required for %s. Install it with install.packages('%s').",
                         pkg, context, pkg), class = "scgeo_missing_pkg")
  }
}
