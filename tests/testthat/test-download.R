listing_html <- '
<a href="../">Parent Directory</a>
<a href="GSE164897_barcodes.tsv.gz">GSE164897_barcodes.tsv.gz</a>
<a href="GSE164897_features.tsv.gz">GSE164897_features.tsv.gz</a>
<a href="GSE164897_matrix.mtx.gz">GSE164897_matrix.mtx.gz</a>
<a href="GSE164897_README.txt">GSE164897_README.txt</a>
'

test_that("scgeo_suppl_files lists files with full FTP urls", {
  testthat::local_mocked_bindings(ftp_get_text = function(url) listing_html)
  res <- scgeo_suppl_files("GSE164897")
  expect_equal(nrow(res), 4)
  expect_true(all(startsWith(res$url,
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164897/suppl/")))
  expect_true("GSE164897_matrix.mtx.gz" %in% res$file)
})

test_that("scgeo_download filters by pattern and writes files", {
  tmp <- withr::local_tempdir()
  testthat::local_mocked_bindings(
    ftp_get_text = function(url) listing_html,
    ftp_download = function(url, path) { writeLines("x", path); path }
  )
  paths <- scgeo_download("GSE164897", pattern = "matrix|barcodes|features", dest = tmp)
  expect_length(paths, 3)
  expect_true(all(file.exists(paths)))
  expect_false(any(grepl("README", paths)))
})

test_that("scgeo_download errors when nothing matches", {
  testthat::local_mocked_bindings(ftp_get_text = function(url) listing_html)
  expect_error(scgeo_download("GSE164897", pattern = "nomatch"),
               class = "scgeo_no_files")
})

test_that("scgeo_download skips re-download when cached", {
  tmp <- withr::local_tempdir()
  calls <- 0
  testthat::local_mocked_bindings(
    ftp_get_text = function(url) listing_html,
    ftp_download = function(url, path) { calls <<- calls + 1; writeLines("x", path); path }
  )
  scgeo_download("GSE164897", pattern = "README", dest = tmp)
  scgeo_download("GSE164897", pattern = "README", dest = tmp)  # cached
  expect_equal(calls, 1)
})

test_that("cache dir helpers work", {
  withr::local_envvar(SCGEO_CACHE = withr::local_tempdir())
  d <- scgeo_cache_dir()
  expect_true(dir.exists(d))
  expect_true(scgeo_clear_cache())
})
