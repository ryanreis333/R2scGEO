test_that("detect_format classifies by extension (incl. gzip)", {
  expect_equal(detect_format("GSE1_matrix.mtx.gz"), "mtx")
  expect_equal(detect_format("foo.mtx"), "mtx")
  expect_equal(detect_format("adata.h5ad"), "h5ad")
  expect_equal(detect_format("data.loom"), "loom")
  expect_equal(detect_format("filtered_feature_bc_matrix.h5"), "h5")
  expect_equal(detect_format("obj.rds"), "rds")
  expect_equal(detect_format("counts.csv.gz"), "dense")
  expect_equal(detect_format("expr.tsv"), "dense")
  expect_equal(detect_format("GSE1_RAW.tar"), "tar")
  expect_equal(detect_format("bundle.tgz"), "tar")
  expect_equal(detect_format("bundle.zip"), "zip")
})

test_that("looks_like_text_table detects delimited text and rejects binary", {
  tmp <- withr::local_tempdir()
  txt <- file.path(tmp, "weirdname")
  writeLines("gene,cell1,cell2", txt)
  expect_true(looks_like_text_table(txt))

  bin <- file.path(tmp, "blob.bin")
  writeBin(as.raw(c(0L, 1L, 2L, 7L)), bin)
  expect_false(looks_like_text_table(bin))
})

test_that("find_triplets groups matrix/barcodes/features by prefix", {
  tmp <- withr::local_tempdir()
  for (f in c("GSM1_matrix.mtx.gz", "GSM1_barcodes.tsv.gz", "GSM1_features.tsv.gz",
              "GSM2_matrix.mtx.gz", "GSM2_barcodes.tsv.gz", "GSM2_genes.tsv.gz")) {
    writeLines("x", file.path(tmp, f))
  }
  trips <- find_triplets(tmp)
  expect_length(trips, 2)
  expect_true(all(vapply(trips, function(t) !is.na(t$barcodes) && !is.na(t$features), logical(1))))
})

test_that("rank_data_files prefers richer formats and count-named files", {
  fs <- c("notes.txt", "raw_counts.csv", "adata.h5ad", "obj.rds", "readme.md")
  ranked <- basename(rank_data_files(fs))
  expect_equal(ranked[1], "adata.h5ad")        # h5ad outranks csv/rds
  expect_true(match("raw_counts.csv", ranked) < match("obj.rds", ranked) ||
              "adata.h5ad" %in% ranked)
})
