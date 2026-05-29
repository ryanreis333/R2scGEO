# ---------------------------------------------------------------------------
# Format detection and file grouping
# ---------------------------------------------------------------------------

# Strip a trailing compression suffix, returning name + whether it was gzipped.
.decompress_name <- function(name) {
  m <- regmatches(name, regexpr("\\.(gz|bz2|xz|zst)$", name, ignore.case = TRUE))
  list(name = sub("\\.(gz|bz2|xz|zst)$", "", name, ignore.case = TRUE),
       compressed = length(m) > 0)
}

#' Identify the storage format of a single file
#'
#' @param path File path.
#' @return One of "mtx", "h5ad", "loom", "h5", "rds", "tar", "zip", "dense",
#'   or "unknown".
#' @keywords internal
detect_format <- function(path) {
  base <- tolower(basename(path))

  # Archives first (may themselves carry .gz).
  if (grepl("\\.(tar\\.gz|tgz|tar)$", base)) return("tar")
  if (grepl("\\.zip$", base)) return("zip")

  inner <- .decompress_name(base)$name   # name without .gz/.bz2/...

  if (grepl("\\.mtx$", inner) || grepl("matrix\\.mtx", base)) return("mtx")
  if (grepl("\\.h5ad$", inner)) return("h5ad")
  if (grepl("\\.loom$", inner)) return("loom")
  if (grepl("\\.(h5|hdf5)$", inner)) return("h5")
  if (grepl("\\.rds$", inner)) return("rds")
  if (grepl("\\.(csv|tsv|txt|tab|counts|expr|matrix)$", inner)) return("dense")

  # Content sniff for extensionless / oddly named text tables.
  if (looks_like_text_table(path)) return("dense")
  "unknown"
}

# Peek at the first line; treat as a dense table if it is delimited text.
looks_like_text_table <- function(path, n = 2) {
  if (!file.exists(path)) return(FALSE)
  con <- tryCatch(
    if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rb") else file(path, "rb"),
    error = function(e) NULL
  )
  if (is.null(con)) return(FALSE)
  on.exit(close(con), add = TRUE)
  bytes <- tryCatch(readBin(con, what = "raw", n = 4096), error = function(e) raw(0))
  if (length(bytes) == 0 || any(bytes == as.raw(0))) return(FALSE)
  text <- rawToChar(bytes)
  lines <- strsplit(text, "\r\n|\n|\r", perl = TRUE)[[1]]
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) return(FALSE)
  head <- utils::head(lines, n)
  if (any(grepl("[\\x01-\\x08\\x0e-\\x1f]", head, perl = TRUE))) return(FALSE)
  grepl("[,\t; ]", head[1])
}

# Group the files in a directory into 10x triplets keyed by shared prefix.
# Returns a list of lists, each with matrix/barcodes/features paths.
find_triplets <- function(dir) {
  fs <- list.files(dir, full.names = TRUE, recursive = TRUE)
  mtx <- grep("matrix\\.mtx(\\.gz)?$", fs, value = TRUE, ignore.case = TRUE)
  if (length(mtx) == 0) return(list())

  lapply(mtx, function(m) {
    prefix <- sub("matrix\\.mtx(\\.gz)?$", "", basename(m), ignore.case = TRUE)
    sib <- fs[dirname(fs) == dirname(m)]
    match_sib <- function(suffixes) {
      for (suf in suffixes) {
        pat <- paste0("^", escape_re(prefix), suf, "(\\.gz)?$")
        hit <- sib[grepl(pat, basename(sib), ignore.case = TRUE)]
        if (length(hit)) return(hit[1])
        # fall back: any file in the folder matching the suffix
        hit <- sib[grepl(paste0(suf, "(\\.gz)?$"), basename(sib), ignore.case = TRUE)]
        if (length(hit)) return(hit[1])
      }
      NA_character_
    }
    list(
      matrix   = m,
      barcodes = match_sib("barcodes\\.tsv"),
      features = match_sib(c("features\\.tsv", "genes\\.tsv"))
    )
  })
}

escape_re <- function(x) gsub("([.\\^$|()\\[\\]{}*+?])", "\\\\\\1", x, perl = TRUE)

# Extract a tar/zip archive into a fresh temp dir and return that dir.
extract_archive <- function(path, fmt = detect_format(path)) {
  out <- tempfile(pattern = paste0("R2scGEO_", gsub("[^A-Za-z0-9]", "_", basename(path)), "_"))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  if (fmt == "zip") {
    utils::unzip(path, exdir = out)
  } else {
    utils::untar(path, exdir = out)
  }
  out
}
