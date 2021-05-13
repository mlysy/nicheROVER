#' Random draws from the posterior distribution with Normal-Inverse-Wishart (NIW) prior.
#'
#' Given iid \eqn{d}-dimensional niche indicators  \eqn{X = (X_1,\ldots,X_N)} with \eqn{X_i \sim N(\mu, \Sigma)}, this function generates random draws from \eqn{p(\mu,\Sigma | X)} for the Normal-Inverse-Wishart (NIW) prior.
#'
#' @details The NIW distribution \eqn{p(\mu, \Sigma | \lambda, \kappa, \Psi, \nu)} is defined as
#' \deqn{
#' \Sigma \sim W^{-1}(\Psi, \nu), \quad \mu | \Sigma \sim N(\lambda, \Sigma/\kappa).
#' }
#' The default value `kappa = 0` uses the Lebesque prior on \eqn{\mu}: \eqn{p(\mu) \propto 1}.
#'
#' The default value `Psi = 0` uses the scale-invariant prior on \eqn{\Sigma}: \eqn{p(\Sigma) \propto |\Sigma|^{-(\nu+d+1)/2}}.
#'
#' The default value `nu = ncol(X)+1` for `kappa = 0` and `Psi = 0` makes \eqn{E[\mu|X]=\code{colMeans(X)}} and \eqn{E[\Sigma | X]=\code{var(X)}}.
#' @param nsamples The number of posterior draws.
#' @param X A data matrix with observations along the rows.
#' @param lambda Location parameter. See 'Details'.
#' @param kappa Scale parameter. Defaults to `kappa = 0`.  See 'Details'.
#' @param Psi Scale matrix. Defaults to `Psi = 0`.  See Details.
#' @param nu Degrees of freedom. Defaults to `nu = ncol(X)+1`.  See Details.
#' @return Returns a list with elements `mu` and `Sigma` of sizes `c(nsamples, length(lambda))` and `c(dim(Psi), nsamples)`.
#' @seealso [rniw()], [niiw.post()].
#' @example examples/niw.post.R
#' @export
niw.post <- function(nsamples, X, lambda, kappa, Psi, nu) {
  par <- niw.coeffs(X, lambda, kappa, Psi, nu)
  rniw(nsamples, par$lambda, par$kappa, par$Psi, par$nu)
}
