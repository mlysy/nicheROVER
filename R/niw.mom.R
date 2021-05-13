#' Mean and variance of the Normal-Inverse-Wishart distribution.
#'
#' This function computes the mean and variance of the Normal-Inverse-Wishart (NIW) distribution.  Can be used to very quickly compute Bayesian point estimates for the conjugate NIW prior.
#'
#' @details The NIW distribution \eqn{p(\mu, \Sigma | \lambda, \kappa, \Psi, \nu)} is defined as
#' \deqn{
#' \Sigma \sim W^{-1}(\Psi, \nu), \quad \mu | \Sigma \sim N(\lambda, \Sigma/\kappa).
#' }
#' Note that cov\eqn{(\mu, \Sigma) = 0}.
#'
#' @param lambda location parameter. See Details.
#' @param kappa scale parameter. See Details.
#' @param Psi scale matrix.  See Details
#' @param nu degrees of freedom.  See Details.
#'
#' @return Returns a list with elements `mu` and `Sigma`, each containing lists with elements `mean` and `var`.  For `mu` these elements are of size `length(lambda)` and `c(length(lambda),length(lambda))`.  For `Sigma` they are of size `dim(Psi)` and `c(dim(Psi), dim(Psi))`, such that cov\eqn{(\Sigma_{ij}, \Sigma_{kl})=}`Sigma$var[i,j,k,l]`.
#'
#' @seealso [rniw()], [niw.coeffs()], [niw.post()].
#' @example examples/niw.mom.R
#' @export
niw.mom <- function(lambda, kappa, Psi, nu) {
  d <- length(lambda)
  b <- nu-d
  mu.mean <- lambda
  Sigma.mean <- Psi/(b-1)
  mu.var <- Sigma.mean/kappa
  Sigma.var <- Psi %o% Psi
  Sigma.var <- 2 * Sigma.var + (b-1) * (aperm(Sigma.var, c(1,3,2,4)) +
                                        aperm(Sigma.var, c(1,4,3,2)))
  Sigma.var <- Sigma.var/(b*(b-1)^2*(b-3))
  list(mu = list(mean = mu.mean, var = mu.var),
       Sigma = list(mean = Sigma.mean, var = Sigma.var))
}
