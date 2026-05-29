#' @keywords internal
"_PACKAGE"

eutils_base <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
geo_ftp_base <- "https://ftp.ncbi.nlm.nih.gov/geo"

# Shared query params NCBI asks all callers to send.
.ncbi_params <- function() {
  p <- list(tool = "R2scGEO", email = Sys.getenv("ENTREZ_EMAIL", ""))
  key <- Sys.getenv("ENTREZ_KEY", "")
  if (nzchar(key)) p$api_key <- key
  p[nzchar(unlist(p))]
}

# Single choke point for E-utilities GET requests. Returns parsed object.
# `parser` is one of "json" or "text". Isolated so tests can stub it.
eutils_get <- function(endpoint, query, parser = "json") {
  req <- httr2::request(eutils_base)
  req <- httr2::req_url_path_append(req, endpoint)
  req <- httr2::req_url_query(req, !!!query, !!!.ncbi_params())
  req <- httr2::req_user_agent(req, "R2scGEO (https://github.com/ryanreis333/R2scGEO)")
  req <- httr2::req_retry(req, max_tries = 3)
  req <- httr2::req_throttle(req, rate = 3)
  resp <- httr2::req_perform(req)
  if (identical(parser, "json")) {
    jsonlite::fromJSON(httr2::resp_body_string(resp), simplifyVector = FALSE)
  } else {
    httr2::resp_body_string(resp)
  }
}

# Validate / normalise a GEO accession like "GSE164897".
parse_accession <- function(acc) {
  acc <- toupper(trimws(acc))
  if (!grepl("^GS[EM][0-9]+$", acc)) {
    rlang::abort(sprintf("'%s' is not a valid GSE/GSM accession.", acc),
                 class = "scgeo_bad_accession")
  }
  acc
}

# Map an accession to its GEO FTP supplementary directory.
# GSE164897 -> .../series/GSE164nnn/GSE164897/suppl/
# GSM5678   -> .../samples/GSM5nnn/GSM5678/suppl/
accession_ftp_dir <- function(acc) {
  acc <- parse_accession(acc)
  prefix <- substr(acc, 1, 3)             # GSE or GSM
  num <- substr(acc, 4, nchar(acc))       # numeric portion as string
  folder <- if (nchar(num) <= 3) {
    paste0(prefix, "nnn")
  } else {
    paste0(prefix, substr(num, 1, nchar(num) - 3), "nnn")
  }
  sub <- if (prefix == "GSE") "series" else "samples"
  sprintf("%s/%s/%s/%s/suppl/", geo_ftp_base, sub, folder, acc)
}

# Build an Entrez `term` from user-facing search arguments.
build_term <- function(query, organism = NULL, single_cell = TRUE) {
  parts <- character(0)
  if (!is.null(query) && nzchar(query)) parts <- c(parts, sprintf("(%s)", query))
  if (single_cell) {
    parts <- c(parts, '("single cell"[All Fields] OR "single-cell"[All Fields] OR scRNA[All Fields])')
  }
  if (!is.null(organism) && nzchar(organism)) {
    parts <- c(parts, sprintf("%s[Organism]", organism))
  }
  # restrict to series-level records
  parts <- c(parts, "GSE[Entry Type]")
  paste(parts, collapse = " AND ")
}
