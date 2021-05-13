#' Random draws from a Normal-Inverse-Wishart distribution.
#'
#'  Generates random draws from a Normal-Inverse-Wishart (NIW) distribution. Can be used to compare prior to posterior parameter distributions.
#'
#' @details The NIW distribution \eqn{p(\mu, \Sigma | \lambda, \kappa, \Psi, \nu)} is defined as
#' \deqn{
#' \Sigma \sim W^{-1}(\Psi, \nu), \quad \mu | \Sigma \sim N(\lambda, \Sigma/\kappa).
#' }
#'
#' @param n number of samples to draw.
#' @param lambda location parameter. See Details.
#' @param kappa scale parameter. See Details.
#' @param Psi scale matrix.  See Details
#' @param nu degrees of freedom.  See Details.
#' @return Returns a list with elements `mu` and `Sigma` of sizes `c(n,length(lambda))` and `c(nrow(Psi),ncol(Psi),n)`.
#'
#' @example examples/rniw.R
#' @seealso [rwish()], [niw.mom()], [niw.coeffs()].
#' @export
rniw <- function(n, lambda, kappa, Psi, nu) {
  d <- length(lambda)
  Sigma <- rwish(n, Psi, nu, inv = TRUE)
  mu <- matrix(NA, n, d)
  colnames(mu) <- names(lambda)
  for(ii in 1:n) {
    mu[ii,] <- .rmvn(mean = lambda, sigma = Sigma[,,ii]/kappa)
  }
  list(mu = mu, Sigma = Sigma)
}
