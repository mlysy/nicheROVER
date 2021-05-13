#' Uniform sampling from an elliptical niche region.
#'
#' @param n Number of random draws.
#' @param mu Mean vector.
#' @param Sigma Variance matrix.
#' @param alpha Probabilistic niche size
#' @return IID draws from a uniform distribution on the elliptical niche region.
#' @example examples/niche.runif.R
#' @seealso [ellipse()] and [niche.size()] for the definition of the elliptical niche region.
#' @export
niche.runif <- function(n, mu, Sigma, alpha = .95) {
  k <- length(mu) # number of dimensions
  # radius draws
  R <- sqrt(qchisq(alpha, df = k)) * rbeta(n, shape1 = k, shape2 = 1)
  Z <- matrix(rnorm(n*k), k, n) # uniform direction
  X <- R * t(Z)/sqrt(colSums(Z*Z)) # circular problem
  t(t(X %*% chol(Sigma)) + mu) # ellipse problem
}
