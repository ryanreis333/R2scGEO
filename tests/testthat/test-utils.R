test_that("parse_accession accepts valid and rejects invalid accessions", {
  expect_equal(parse_accession(" gse164897 "), "GSE164897")
  expect_equal(parse_accession("GSM5678"), "GSM5678")
  expect_error(parse_accession("GSE"), class = "scgeo_bad_accession")
  expect_error(parse_accession("ABC123"), class = "scgeo_bad_accession")
})

test_that("accession_ftp_dir builds correct GEO FTP paths", {
  expect_equal(
    accession_ftp_dir("GSE164897"),
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164897/suppl/"
  )
  expect_equal(
    accession_ftp_dir("GSE123"),
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE123/suppl/"
  )
  expect_equal(
    accession_ftp_dir("GSE1234"),
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1234/suppl/"
  )
  expect_equal(
    accession_ftp_dir("GSM5678"),
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5nnn/GSM5678/suppl/"
  )
})

test_that("build_term composes Entrez query parts", {
  t <- build_term("glioma", organism = "Homo sapiens", single_cell = TRUE)
  expect_match(t, "\\(glioma\\)")
  expect_match(t, "Homo sapiens\\[Organism\\]")
  expect_match(t, "single-cell", fixed = TRUE)
  expect_match(t, "GSE\\[Entry Type\\]")

  t2 <- build_term("glioma", single_cell = FALSE)
  expect_false(grepl("single-cell", t2, fixed = TRUE))
})

test_that("parse_apache_listing extracts file names and drops nav links", {
  html <- '<a href="/geo/series/">Parent</a>
           <a href="../">..</a>
           <a href="GSE164897_RAW.tar">GSE164897_RAW.tar</a>
           <a href="GSE164897_counts.mtx.gz">counts</a>
           <a href="https://www.hhs.gov/vulnerability-disclosure-policy/index.html">HHS</a>
           <a href="subdir/">subdir/</a>'
  files <- parse_apache_listing(html)
  expect_setequal(files, c("GSE164897_RAW.tar", "GSE164897_counts.mtx.gz"))
})
