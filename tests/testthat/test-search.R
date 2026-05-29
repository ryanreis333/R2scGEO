# Canned E-utilities responses (shape matches jsonlite::fromJSON(simplifyVector = FALSE)).
fake_esearch <- function(n = 1) {
  list(esearchresult = list(idlist = as.list(rep("200164897", n))))
}
fake_esummary <- function() {
  list(result = list(
    uids = list("200164897"),
    `200164897` = list(
      accession = "GSE164897",
      title = "scRNA-seq of human glioblastoma",
      taxon = "Homo sapiens",
      n_samples = "12",
      gdsType = "Expression profiling by high throughput sequencing",
      PDAT = "2021/03/15",
      summary = "Single-cell transcriptomes from tumor samples."
    )
  ))
}

# Dispatcher that mimics eutils_get based on endpoint.
mock_eutils <- function(endpoint, query, parser = "json") {
  if (grepl("esearch", endpoint)) fake_esearch() else fake_esummary()
}

test_that("scgeo_search returns a tidy one-row tibble", {
  testthat::local_mocked_bindings(eutils_get = mock_eutils)
  res <- scgeo_search("glioblastoma", organism = "Homo sapiens")
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 1)
  expect_equal(res$accession, "GSE164897")
  expect_equal(res$organism, "Homo sapiens")
  expect_equal(res$n_samples, 12L)
  expect_type(res$n_samples, "integer")
})

test_that("scgeo_search returns empty tibble when no ids found", {
  testthat::local_mocked_bindings(
    eutils_get = function(endpoint, query, parser = "json") {
      list(esearchresult = list(idlist = list()))
    }
  )
  res <- scgeo_search("nonexistent")
  expect_equal(nrow(res), 0)
  expect_named(res, c("accession", "title", "organism", "n_samples",
                      "gds_type", "pubdate", "summary"))
})

test_that("scgeo_metadata resolves accession to metadata row", {
  testthat::local_mocked_bindings(eutils_get = mock_eutils)
  res <- scgeo_metadata("GSE164897")
  expect_equal(nrow(res), 1)
  expect_equal(res$title, "scRNA-seq of human glioblastoma")
})

test_that("scgeo_metadata rejects bad accessions", {
  expect_error(scgeo_metadata("not-an-accession"), class = "scgeo_bad_accession")
})
