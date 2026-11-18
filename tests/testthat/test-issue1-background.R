# Regression tests for GH issue #1: ultragsea() failed with a single gene
# set because the background was the intersection of fc and rownames(G)
# (frequently just the gene set itself), collapsing G to a vector and/or a
# zero-variance column. The fix uses the full set of fc genes as background.

test_that("ultragsea works with a single gene set (no crash, finite result)", {
  n_genes <- 1e3L
  set.seed(0L)
  fc <- stats::rnorm(n_genes)
  genes <- paste0("gene", seq_len(n_genes))
  names(fc) <- genes

  gene_sets <- list(set1 = genes[seq_len(10L)])
  G <- gmt2mat(gene_sets)

  res <- ultragsea(G, fc)

  expect_equal(nrow(res), 1)
  expect_false(is.na(res$score))
  expect_false(is.na(res$pval))
  expect_true(is.finite(res$score))
})

test_that("ultragsea works when one gene set is nested in another", {
  n_genes <- 1e3L
  set.seed(0L)
  fc <- stats::rnorm(n_genes)
  genes <- paste0("gene", seq_len(n_genes))
  names(fc) <- genes

  gene_sets <- list(set1 = genes[seq_len(10L)], set2 = genes[seq_len(5L)])
  G <- gmt2mat(gene_sets)

  res <- ultragsea(G, fc)

  expect_equal(nrow(res), 2)
  expect_true(all(is.finite(res$score)))
  expect_true(all(is.finite(res$pval)))
})

test_that("ultragsea returns NA (not NaN/Inf/error) when background has zero variance", {
  # background == the gene set itself -> the membership column is constant,
  # so the correlation is genuinely undefined (0/0), and must surface as NA.
  genes <- paste0("gene", 1:10)
  set.seed(1L)
  fc <- stats::setNames(stats::rnorm(10), genes)
  G <- gmt2mat(list(set1 = genes))

  res <- ultragsea(G, fc)

  expect_true(is.na(res$score))
  expect_true(is.na(res$pval))
  expect_false(is.nan(res$score))
})

test_that("gset.cortest aligns gset rows to FC order even without padding", {
  # regression: sparsematrix_match() must reindex rows/cols even when no
  # gene is missing, otherwise a differently-ordered gset silently
  # misaligns with FC.
  fc <- stats::setNames(c(1, 2, 3, 4, 5), c("g1", "g2", "g3", "g4", "g5"))
  gset <- Matrix::Matrix(c(1, 0, 1, 0, 1), ncol = 1, sparse = TRUE)
  rownames(gset) <- c("g5", "g4", "g3", "g2", "g1")
  colnames(gset) <- "set1"

  res <- gset.cortest(gset, fc)

  expect_equal(as.numeric(res$rho[, 1]), stats::cor(c(1, 0, 1, 0, 1), c(5, 4, 3, 2, 1)))
})
