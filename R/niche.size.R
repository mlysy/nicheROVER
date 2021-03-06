#' Calculate the size of an elliptical niche region.
#'
#' @details For a given niche region \eqn{N_R}, the niche size is defined as the hypervolume of this region: \eqn{N_S = \int_{x \in N_R} d x}.
#'
#' @param Sigma Variance matrix for normally distributed niche axes.
#' @param alpha Probabilistic niche size.
#' @return Hypervolume niche size.
#'
#' @example examples/niche.size.R
#'
#' @export
niche.size <- function(Sigma, alpha = .95) {
  Sigma <- as.matrix(Sigma)
  n <- nrow(Sigma)
  sz <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  sz <- .5 * (sz + n * (log(pi) + log(qchisq(alpha, df = n)))) - lgamma(.5*n+1)
  exp(sz)
}
