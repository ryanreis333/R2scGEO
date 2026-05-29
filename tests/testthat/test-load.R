# Build a tiny real 10x triplet on disk and load it back.
write_triplet <- function(dir, gz = FALSE) {
  ext <- if (gz) ".gz" else ""
  op <- function(f) if (gz) gzfile(f) else file(f)

  mtx <- c(
    "%%MatrixMarket matrix coordinate integer general",
    "%",
    "3 2 3",      # 3 genes, 2 cells, 3 nonzero
    "1 1 5",
    "2 2 3",
    "3 1 1"
  )
  con <- op(file.path(dir, paste0("matrix.mtx", ext))); writeLines(mtx, con); close(con)
  con <- op(file.path(dir, paste0("barcodes.tsv", ext)))
  writeLines(c("AAAC-1", "AAAG-1"), con); close(con)
  con <- op(file.path(dir, paste0("features.tsv", ext)))
  writeLines(c("ENSG1\tGeneA\tGene Expression",
               "ENSG2\tGeneB\tGene Expression",
               "ENSG3\tGeneC\tGene Expression"), con); close(con)
}

test_that("scgeo_load_10x reads a plain triplet into a sparse matrix", {
  tmp <- withr::local_tempdir()
  write_triplet(tmp, gz = FALSE)
  m <- scgeo_load_10x(tmp, as = "matrix")
  expect_s4_class(m, "dgCMatrix")
  expect_equal(dim(m), c(3, 2))
  expect_equal(as.numeric(m["GeneA", "AAAC-1"]), 5)
  expect_equal(rownames(m), c("GeneA", "GeneB", "GeneC"))
  expect_equal(colnames(m), c("AAAC-1", "AAAG-1"))
})

test_that("scgeo_load_10x reads gzipped triplet", {
  tmp <- withr::local_tempdir()
  write_triplet(tmp, gz = TRUE)
  m <- scgeo_load_10x(tmp, as = "matrix")
  expect_equal(dim(m), c(3, 2))
})

test_that("scgeo_load_10x errors without a matrix file", {
  tmp <- withr::local_tempdir()
  expect_error(scgeo_load_10x(tmp, as = "matrix"), class = "scgeo_no_matrix")
})

test_that("seurat output is skipped when Seurat absent", {
  skip_if_not_installed("Seurat")
  tmp <- withr::local_tempdir()
  write_triplet(tmp)
  obj <- scgeo_load_10x(tmp, as = "seurat")
  expect_s4_class(obj, "Seurat")
})
