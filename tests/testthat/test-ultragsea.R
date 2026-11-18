# Unit tests for ultragsea main function

test_that("ultragsea basic functionality works", {
  # Create test data
  set.seed(123)
  fc <- rnorm(100, mean = 0, sd = 1)
  names(fc) <- paste0("gene_", 1:100)

  # Create gene set matrix (binary, sparse)
  G <- matrix(rbinom(100 * 10, 1, 0.3), nrow = 100, ncol = 10)
  rownames(G) <- names(fc)
  colnames(G) <- paste0("geneset_", 1:10)
  G <- Matrix::Matrix(G, sparse = TRUE)

  # Test basic call
  result <- ultragsea(G, fc)

  # Check output structure
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 10)
  expect_equal(ncol(result), 5)
  expect_named(result, c("pathway", "pval", "padj", "score", "size"))
  # Check column types
  expect_type(result$pval, "double")
  expect_type(result$padj, "double")
  expect_type(result$score, "double")
  expect_type(result$size, "integer")
})

test_that("ultragsea with format='as.gsea' returns data.table", {
  # Create test data
  set.seed(456)
  fc <- rnorm(30, mean = 0, sd = 1)
  names(fc) <- paste0("gene_", 1:30)

  # Create gene set matrix
  G <- matrix(rbinom(30 * 3, 1, 0.5), nrow = 30, ncol = 3)
  rownames(G) <- names(fc)
  colnames(G) <- paste0("gs", 1:3)
  G <- Matrix::Matrix(G, sparse = TRUE)

  # Test with as.gsea format
  result <- ultragsea(G, fc, format = "as.gsea")

  # Check output structure
  expect_s3_class(result, "data.table")
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 8)
  expect_named(result,
    c("pathway", "pval", "padj", "log2err", "ES", "NES",
      "size", "leadingEdge"))
})

test_that("ultragsea with different methods", {
  # Create test data
  set.seed(789)
  fc <- rnorm(40, mean = 0, sd = 1)
  names(fc) <- paste0("gene_", 1:40)

  # Create gene set matrix
  G <- matrix(rbinom(40 * 4, 1, 0.35), nrow = 40, ncol = 4)
  rownames(G) <- names(fc)
  colnames(G) <- paste0("set", 1:4)
  G <- Matrix::Matrix(G, sparse = TRUE)

  # Test ztest method
  result_ztest <- ultragsea(G, fc, method = "ztest")
  expect_s3_class(result_ztest, "data.frame")
  expect_equal(nrow(result_ztest), 4)

  # Test ttest method
  result_ttest <- ultragsea(G, fc, method = "ttest")
  expect_s3_class(result_ttest, "data.frame")
  expect_equal(nrow(result_ttest), 4)

  # Test cor method
  result_cor <- ultragsea(G, fc, method = "cor")
  expect_s3_class(result_cor, "data.frame")
  expect_equal(nrow(result_cor), 4)
})

test_that("ultragsea with alpha parameter", {
  # Create test data
  set.seed(202)
  fc <- rnorm(35, mean = 0, sd = 1)
  names(fc) <- paste0("gene_", 1:35)

  # Create gene set matrix
  G <- matrix(rbinom(35 * 3, 1, 0.4), nrow = 35, ncol = 3)
  rownames(G) <- names(fc)
  colnames(G) <- paste0("gs", 1:3)
  G <- Matrix::Matrix(G, sparse = TRUE)

  # Test with different alpha values (ztest: alpha drives Bayesian shrinkage;
  # the default cor method ignores alpha)
  result_low <- ultragsea(G, fc, alpha = 0.1, method = "ztest")
  result_high <- ultragsea(G, fc, alpha = 0.9, method = "ztest")

  # Different alpha values should give different scores
  expect_false(isTRUE(all.equal(result_low$score, result_high$score)))
})

test_that("ultragsea with minLE parameter", {
  # Deterministic data: half of the genes strongly up, half strongly down,
  # so minLE filtering is well-defined without RNG
  fc <- c(rep(5, 13), rep(-5, 12))
  names(fc) <- paste0("g", 1:25)

  # set1 spans both signs (positive correlation with fc); set2 is all-up
  G <- cbind(set1 = c(rep(1, 18), rep(0, 7)),
             set2 = c(rep(1, 13), rep(0, 12)))
  rownames(G) <- names(fc)
  G <- Matrix::Matrix(G, sparse = TRUE)

  result <- ultragsea(G, fc, minLE = 2, format = "as.gsea")

  expect_type(result$leadingEdge, "list")
  expect_equal(length(result$leadingEdge), 2)
  # minLE = 2 keeps only the strongly up-regulated genes in the leading
  # edge; the down-regulated members of set1 are filtered out
  expect_equal(length(result$leadingEdge[[1]]), 13)
  expect_equal(length(result$leadingEdge[[2]]), 13)
})

test_that("ultragsea with gene intersection", {
  # Create test data
  set.seed(321)
  fc <- rnorm(100, mean = 0, sd = 1)
  names(fc) <- paste0("gene_", 1:100)

  # Create gene set matrix with subset of genes
  G <- matrix(rbinom(50 * 5, 1, 0.4), nrow = 50, ncol = 5)
  rownames(G) <- paste0("gene_", 51:100)  # Partial overlap with fc
  colnames(G) <- paste0("geneset_", 1:5)
  G <- Matrix::Matrix(G, sparse = TRUE)

  # Test should handle gene intersection
  result <- ultragsea(G, fc)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 5)
  expect_true(all(result$size <= 50))  # Limited by intersection
})

test_that("ultragsea returns valid p-values", {
  # Create test data
  set.seed(654)
  fc <- rnorm(60, mean = 0, sd = 1)
  names(fc) <- paste0("gene_", 1:60)

  # Create gene set matrix
  G <- matrix(rbinom(60 * 6, 1, 0.3), nrow = 60, ncol = 6)
  rownames(G) <- names(fc)
  colnames(G) <- paste0("set", 1:6)
  G <- Matrix::Matrix(G, sparse = TRUE)

  # Test that p-values are valid
  result <- ultragsea(G, fc)

  expect_true(all(result$pval >= 0 & result$pval <= 1, na.rm = TRUE))
  expect_true(all(result$padj >= 0 & result$padj <= 1, na.rm = TRUE))
  expect_true(all(result$padj >= result$pval - 1e-10, na.rm = TRUE))
})

test_that("ultragsea with invalid method throws error", {
  # Create minimal test data
  fc <- rnorm(20, mean = 0, sd = 1)
  names(fc) <- paste0("gene_", 1:20)

  G <- matrix(rbinom(20 * 2, 1, 0.3), nrow = 20, ncol = 2)
  rownames(G) <- names(fc)
  colnames(G) <- paste0("set", 1:2)

  expect_error(
    ultragsea(G, fc, method = "invalid"),
    "unsupported method"
  )
})

test_that("ultragsea with sparse matrix input", {
  # Create test data
  set.seed(987)
  fc <- rnorm(80, mean = 0, sd = 1)
  names(fc) <- paste0("gene_", 1:80)

  # Create gene set matrix
  G_dense <- matrix(rbinom(80 * 8, 1, 0.25), nrow = 80, ncol = 8)
  rownames(G_dense) <- names(fc)
  colnames(G_dense) <- paste0("set", 1:8)

  # Test sparse matrix input
  G_sparse <- Matrix::Matrix(G_dense, sparse = TRUE)
  result_sparse <- ultragsea(G_sparse, fc, method = "ztest")

  # Test dense matrix input (ztest accepts any matrix)
  result_dense <- ultragsea(G_dense, fc, method = "ztest")

  # Results should be identical
  expect_equal(result_sparse, result_dense)
})