#' Location of the on-disk download cache
#' @return Path to the cache directory (created on first use).
#' @export
scgeo_cache_dir <- function() {
  d <- Sys.getenv("SCGEO_CACHE", tools::R_user_dir("R2scGEO", "cache"))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  d
}

#' Clear the scgeo download cache
#' @param accession Optional accession to clear; if `NULL`, clears everything.
#' @return Invisibly, `TRUE`.
#' @export
scgeo_clear_cache <- function(accession = NULL) {
  root <- scgeo_cache_dir()
  target <- if (is.null(accession)) root else file.path(root, parse_accession(accession))
  if (dir.exists(target)) unlink(target, recursive = TRUE)
  if (is.null(accession) && !dir.exists(root)) dir.create(root, recursive = TRUE)
  invisible(TRUE)
}

#' List supplementary files available for a GEO accession
#'
#' Parses the Apache directory listing of the accession's GEO FTP `suppl/`
#' folder.
#'
#' @param accession A GSE/GSM accession, e.g. "GSE164897".
#' @return A tibble with `file` and `url` columns.
#' @export
scgeo_suppl_files <- function(accession) {
  dir_url <- accession_ftp_dir(accession)
  html <- ftp_get_text(dir_url)
  files <- parse_apache_listing(html)
  tibble::tibble(file = files, url = paste0(dir_url, files))
}

#' Download supplementary files for a GEO accession
#'
#' @param accession A GSE/GSM accession.
#' @param pattern Optional regex; only files whose name matches are downloaded.
#' @param dest Destination directory (defaults to a per-accession cache folder).
#' @param overwrite Re-download even if a cached copy exists (default `FALSE`).
#' @return A character vector of local file paths.
#' @export
scgeo_download <- function(accession, pattern = NULL, dest = NULL, overwrite = FALSE) {
  accession <- parse_accession(accession)
  files <- scgeo_suppl_files(accession)
  if (!is.null(pattern)) files <- files[grepl(pattern, files$file), , drop = FALSE]
  if (nrow(files) == 0) {
    rlang::abort("No supplementary files matched.", class = "scgeo_no_files")
  }
  if (is.null(dest)) dest <- file.path(scgeo_cache_dir(), accession)
  if (!dir.exists(dest)) dir.create(dest, recursive = TRUE)

  vapply(seq_len(nrow(files)), function(i) {
    local <- file.path(dest, files$file[i])
    if (overwrite || !file.exists(local)) ftp_download(files$url[i], local)
    local
  }, character(1))
}

# --- internal helpers -------------------------------------------------------

ftp_get_text <- function(url) {
  req <- httr2::request(url)
  req <- httr2::req_user_agent(req, "R2scGEO")
  req <- httr2::req_retry(req, max_tries = 3)
  httr2::resp_body_string(httr2::req_perform(req))
}

ftp_download <- function(url, path) {
  req <- httr2::request(url)
  req <- httr2::req_user_agent(req, "R2scGEO")
  req <- httr2::req_retry(req, max_tries = 3)
  httr2::req_perform(req, path = path)
  path
}

# Extract file names from an Apache/NCBI directory listing.
parse_apache_listing <- function(html) {
  hrefs <- regmatches(html, gregexpr('href="([^"?][^"]*)"', html, perl = TRUE))[[1]]
  hrefs <- sub('^href="', "", sub('"$', "", hrefs))
  # Keep only files linked relative to the current suppl/ folder.
  hrefs <- hrefs[!grepl("^/|^\\.\\.|/$|^[A-Za-z][A-Za-z0-9+.-]*:", hrefs)]
  unique(hrefs)
}
