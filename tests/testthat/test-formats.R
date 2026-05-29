make_triplet <- function(dir, prefix = "", gz = TRUE) {
  ext <- if (gz) ".gz" else ""
  op <- function(f) if (gz) gzfile(f) else file(f)
  w <- function(name, lines) { con <- op(file.path(dir, paste0(prefix, name, ext)))
                               writeLines(lines, con); close(con) }
  w("matrix.mtx", c("%%MatrixMarket matrix coordinate integer general", "%",
                    "3 2 3", "1 1 5", "2 2 3", "3 1 1"))
  w("barcodes.tsv", c("AAAC-1", "AAAG-1"))
  w("features.tsv", c("ENSG1\tGeneA\tGene Expression",
                      "ENSG2\tGeneB\tGene Expression",
                      "ENSG3\tGeneC\tGene Expression"))
}

test_that("scgeo_load reads a dense csv (genes x cells)", {
  tmp <- withr::local_tempdir()
  f <- file.path(tmp, "counts.csv")
  writeLines(c("gene,c1,c2,c3", "GeneA,5,0,2", "GeneB,0,3,0"), f)
  m <- scgeo_load(f, as = "matrix")
  expect_s4_class(m, "dgCMatrix")
  expect_equal(dim(m), c(2, 3))
  expect_equal(rownames(m), c("GeneA", "GeneB"))
  expect_equal(as.numeric(m["GeneA", "c1"]), 5)
})

test_that("scgeo_load reads a gzipped dense table", {
  tmp <- withr::local_tempdir()
  f <- file.path(tmp, "expr.tsv.gz")
  con <- gzfile(f); writeLines(c("gene\tc1\tc2", "GeneA\t1\t2", "GeneB\t3\t4"), con); close(con)
  m <- scgeo_load(f)
  expect_equal(dim(m), c(2, 2))
})

test_that("scgeo_load reads an .rds matrix and passes through Seurat objects", {
  tmp <- withr::local_tempdir()
  mat <- Matrix::Matrix(c(1, 0, 2, 3), 2, 2, sparse = TRUE)
  rownames(mat) <- c("GeneA", "GeneB"); colnames(mat) <- c("c1", "c2")
  f <- file.path(tmp, "obj.rds"); saveRDS(mat, f)
  m <- scgeo_load(f, as = "matrix")
  expect_s4_class(m, "dgCMatrix")
  expect_equal(dim(m), c(2, 2))
})

test_that("scgeo_load auto-finds a triplet in a directory", {
  tmp <- withr::local_tempdir()
  make_triplet(tmp)
  m <- scgeo_load(tmp, as = "matrix")
  expect_equal(dim(m), c(3, 2))
})

test_that("scgeo_load extracts a tar archive and reads contents", {
  tmp <- withr::local_tempdir()
  inner <- file.path(tmp, "data"); dir.create(inner)
  make_triplet(inner, gz = FALSE)
  tarf <- file.path(tmp, "GSE1_RAW.tar")
  withr::with_dir(inner, utils::tar(tarf, files = list.files(inner)))
  m <- scgeo_load(tarf, as = "matrix")
  expect_equal(dim(m), c(3, 2))
})

test_that("multi-sample directory returns a named list", {
  tmp <- withr::local_tempdir()
  make_triplet(tmp, prefix = "GSM1_")
  make_triplet(tmp, prefix = "GSM2_")
  res <- scgeo_load(tmp, as = "matrix")
  expect_type(res, "list")
  expect_length(res, 2)
  expect_true(all(vapply(res, function(x) methods::is(x, "dgCMatrix"), logical(1))))
})

test_that("sample selector picks one triplet from many", {
  tmp <- withr::local_tempdir()
  make_triplet(tmp, prefix = "GSM1_")
  make_triplet(tmp, prefix = "GSM2_")
  m <- scgeo_load(tmp, as = "matrix", sample = "GSM2")
  expect_s4_class(m, "dgCMatrix")
})

test_that("unknown format raises a clear error", {
  tmp <- withr::local_tempdir()
  f <- file.path(tmp, "mystery.bin")
  writeBin(as.raw(0:9), f)
  expect_error(scgeo_load(f), class = "scgeo_unknown_format")
})

# --- optional formats (need extra packages) --------------------------------

test_that("10x .h5 loads when Seurat or hdf5r present", {
  skip_if_not(requireNamespace("Seurat", quietly = TRUE) ||
              requireNamespace("hdf5r", quietly = TRUE))
  skip("requires a real CellRanger .h5 fixture; covered by integration tests")
})

test_that(".h5ad loads when a backend is present", {
  skip_if_not(requireNamespace("zellkonverter", quietly = TRUE) ||
              requireNamespace("anndata", quietly = TRUE) ||
              requireNamespace("hdf5r", quietly = TRUE))
  skip("requires a real .h5ad fixture; covered by integration tests")
})
