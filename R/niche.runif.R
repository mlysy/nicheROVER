#' Uniform sampling from elliptical niche region.
#'
#' @param n number of random draws
#' @param mu mean vector
#' @param Sigma variance matrix
#' @param alpha probabilistic niche size
#' @return iid draws from uniform distribution on niche region.
#' @examples
#' # 2d example
#' V <- crossprod(matrix(rnorm(4),2,2))
#' mu <- rnorm(2)
#' plot(ellipse(mu, V), type = "l")
#' points(niche.runif(1e4, mu, V), col = "brown", pch = ".")
#' @export
niche.runif <- function(n, mu, Sigma, alpha = .95) {
  k <- length(mu) # number of dimensions
  # radius draws
  R <- sqrt(qchisq(alpha, df = k)) * rbeta(n, shape1 = k, shape2 = 1)
  Z <- matrix(rnorm(n*k), k, n) # uniform direction
  X <- R * t(Z)/sqrt(colSums(Z*Z)) # circular problem
  t(t(X %*% chol(Sigma)) + mu) # ellipse problem
}
