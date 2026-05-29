#' Search GEO for single-cell RNA-seq series
#'
#' Runs an Entrez `esearch` against the GEO DataSets (`gds`) database and
#' returns a tidy table of matching series with their key metadata.
#'
#' @param query Free-text query, e.g. "glioblastoma" or "T cell exhaustion".
#' @param organism Optional organism filter, e.g. "Homo sapiens".
#' @param max_results Maximum number of series to return (default 25).
#' @param single_cell If `TRUE` (default), restrict results to single-cell studies.
#'
#' @return A [tibble][tibble::tibble] with one row per series: `accession`,
#'   `title`, `organism`, `n_samples`, `gds_type`, `pubdate`, `summary`.
#' @export
#' @examples
#' \dontrun{
#' scgeo_search("glioblastoma", organism = "Homo sapiens", max_results = 10)
#' }
scgeo_search <- function(query, organism = NULL, max_results = 25, single_cell = TRUE) {
  term <- build_term(query, organism, single_cell)
  es <- eutils_get("esearch.fcgi",
                   list(db = "gds", term = term, retmax = max_results, retmode = "json"))
  ids <- unlist(es$esearchresult$idlist)
  if (length(ids) == 0) {
    return(empty_result())
  }
  summary_to_tibble(fetch_summaries(ids))
}

#' Fetch metadata for one or more GEO accessions
#'
#' @param accession Character vector of GSE/GSM accessions, e.g. "GSE164897".
#' @return A tibble, one row per accession (same columns as [scgeo_search()]).
#' @export
scgeo_metadata <- function(accession) {
  accession <- vapply(accession, parse_accession, character(1))
  # esearch each accession to obtain its internal UID, then esummary.
  ids <- vapply(accession, function(a) {
    es <- eutils_get("esearch.fcgi",
                     list(db = "gds", term = sprintf("%s[Accession]", a), retmode = "json"))
    id <- unlist(es$esearchresult$idlist)
    if (length(id) == 0) NA_character_ else id[1]
  }, character(1))
  ids <- ids[!is.na(ids)]
  if (length(ids) == 0) return(empty_result())
  summary_to_tibble(fetch_summaries(unname(ids)))
}

# --- internal helpers -------------------------------------------------------

fetch_summaries <- function(ids) {
  res <- eutils_get("esummary.fcgi",
                    list(db = "gds", id = paste(ids, collapse = ","), retmode = "json"))
  docs <- res$result
  docs[setdiff(names(docs), "uids")]
}

summary_to_tibble <- function(docs) {
  if (length(docs) == 0) return(empty_result())
  rows <- lapply(docs, function(d) {
    tibble::tibble(
      accession = chr(d$accession),
      title     = chr(d$title),
      organism  = chr(d$taxon),
      n_samples = int(d$n_samples),
      gds_type  = chr(d$gdsType),
      pubdate   = chr(d$PDAT),
      summary   = chr(d$summary)
    )
  })
  do.call(rbind, rows)
}

empty_result <- function() {
  tibble::tibble(
    accession = character(0), title = character(0), organism = character(0),
    n_samples = integer(0), gds_type = character(0), pubdate = character(0),
    summary = character(0)
  )
}

chr <- function(x) if (is.null(x) || length(x) == 0) NA_character_ else as.character(x)[1]
int <- function(x) {
  if (is.null(x) || length(x) == 0 || !nzchar(as.character(x)[1])) return(NA_integer_)
  suppressWarnings(as.integer(x[1]))
}
