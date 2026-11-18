
#' Calculate gene set correlation
#'
#' Compute correlation between a foldchange vector/matrix and gene sets
#'
#' @param FC Numeric vector or matrix of (log)foldchanges, with genes as row names
#' @param gset Numeric matrix of gene sets, with genes as row/column names
#' @param compute.p Logical indicating whether to compute p-values
#' @param use.rank Logical indicating whether to rank transform FC before correlation
#'
#' @return Named list with components:
#' \itemize{
#'  \item rho - Matrix of correlation coefficients between FC and gset
#'  \item p.value - Matrix of p-values for correlation (if compute.p = TRUE)
#'  \item q.value - Matrix of FDR adjusted p-values (if compute.p = TRUE)
#' }
#'
#' @details This function calculates sparse rank correlation between FC and each
#' column of gset using \code{qlcMatrix::corSparse()}. It handles missing values in
#' FC by computing column-wise correlations.
#'
#' P-values are computed from statistical distribution
#'
#' @examples
#' \dontrun{
#' library(playbase)
#' ranks <- sample(1:10000, 1000, replace = TRUE)
#' names(ranks) <- replicate(1000, paste(sample(LETTERS, 4, replace = TRUE), collapse = ""))
#' genesets <- matrix(rnorm(1000 * 20), ncol = 20)
#' rownames(genesets) <- names(ranks)
#'
#' gset.rankcor(ranks, genesets, compute.p = TRUE)
#' }
#' @export
gset.cortest <- function(gset, FC, compute.p = FALSE, use.rank = FALSE,
                         corshrink = 0) {
  if (ncol(gset) == 0 || NCOL(FC) == 0) {
    if (ncol(gset) == 0) stop("gset has zero columns")
    if (NCOL(FC) == 0) stop("FC has zero columns")
  }

  if (!inherits(gset, "Matrix")) stop("gset must be a matrix")

  is.vec <- (NCOL(FC) == 1 && !inherits(FC, c("matrix", "Matrix")))
  if (is.vec && is.null(names(FC))) stop("rank vector must be named")
  if (!is.vec && is.null(rownames(FC))) stop("rank matrix must have rownames")
  if (is.vec) FC <- matrix(FC, ncol = 1, dimnames = list(names(FC), "FC"))
  n1 <- sum(rownames(FC) %in% colnames(gset), na.rm = TRUE)
  n2 <- sum(rownames(FC) %in% rownames(gset), na.rm = TRUE)
  if (n1 > n2) gset <- Matrix::t(gset)

  # align gset to the background given by FC (see ultragsea(), #1)
  gset <- sparsematrix_match(gset, rows=rownames(FC), cols=NULL)


  if (use.rank) {
    if (inherits(FC, "dgCMatrix")) {
      ## for sparse dgCMatrix
      FC <- sparseMatrixStats::colRanks(FC, na.last = "keep", ties.method = "random", preserveShape = TRUE)
    } else {
      FC <- matrixStats::colRanks(FC, na.last = "keep", ties.method = "random", preserveShape = TRUE)
    }
  }

  ## two cases: (1) in case no missing values, just use corSparse on
  ## whole matrix. (2) in case the FC matrix has missing values, we
  ## must proceed 1-column at time and do reduced corSparse on
  ## intersection of genes.
  rho1 <- cor_sparse_matrix(gset, FC)

  rownames(rho1) <- colnames(gset)
  colnames(rho1) <- colnames(FC)
  rho1[is.nan(rho1)] <- NA ## ??

  if(corshrink > 0) {
    ## Shrinkage to penalize small gene sets.
    gsize <- Matrix::colSums(gset!=0)
    rho1 <- rho1 / (1 + (corshrink / gsize)) 
  }
  
  ## compute p-value
  if (compute.p) {
    pv <- apply(rho1, 2, function(x) cor_pvalue(x, n = nrow(FC)))
    pv[is.nan(pv)] <- NA ## ??
    pv <- matrix(pv, nrow = nrow(rho1), ncol = ncol(rho1),
      dimnames = dimnames(rho1))   ## keep dims for single set
    qv <- apply(pv, 2, stats::p.adjust, method = "fdr")
    qv <- matrix(qv, nrow = nrow(pv), ncol = ncol(pv), dimnames = dimnames(pv))
    df <- list(rho = rho1, p.value = pv, q.value = qv)
  } else {
    df <- list(rho = rho1, p.value = NA, q.value = NA)
  }
  df
}

#' Calculate p-value from correlation rho 
#'
#' @param x  correlation value
#' @param n  degrees of freendom
#' @return p valu
#'
cor_pvalue <- function(x, n) 2 * stats::pnorm(-abs(x / ((1 - x**2) / (n - 2))**0.5))

#' Calculate sparse correlation matrix handling missing values
#'
#' @param G Sparse matrix containing gene sets
#' @param mat Matrix of values
#' @return Correlation matrix between G and mat
#' @details If mat has no missing values, calculates correlation directly using corSparse.
#' Otherwise computes column-wise correlations only using non-missing values.
#'
cor_sparse_matrix <- function(G, mat) {
  if (sum(is.na(mat)) == 0) {
    cor_matrix <- qlcMatrix::corSparse(G, mat)
  } else {
    message("matrix has missing values: computing column-wise reduced cor")
    corSparse.vec <- function(X, y) {
      jj <- which(!is.na(y))
      qlcMatrix::corSparse(X[jj, , drop = FALSE], cbind(y[jj]))
    }
    cor_matrix <- lapply(1:ncol(mat), function(i) corSparse.vec(G, mat[, i]))
    cor_matrix <- do.call(cbind, cor_matrix)
  }
  return(cor_matrix)
}

