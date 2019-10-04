#' Random draws from the posterior distribution with Normal-Inverse-Wishart (NIW) prior.
#'
#' Given iid \eqn{d}-dimensional niche indicators  \eqn{X = (X_1,\ldots,X_N)} with \eqn{X_i \sim N(\mu, \Sigma)}, this function generates random draws from \eqn{p(\mu,\Sigma | X)} for the Normal-Inverse-Wishart (NIW) prior.
#'
#' @details The NIW distribution \eqn{p(\mu, \Sigma | \lambda, \kappa, \Psi, \nu)} is defined as
#' \deqn{
#' \Sigma \sim W^{-1}(\Psi, \nu), \quad \mu | \Sigma \sim N(\lambda, \Sigma/\kappa).
#' }
#' The default value \code{kappa = 0} uses the Lebesque prior on \eqn{\mu}: \eqn{p(\mu) \propto 1}.
#'
#' The default value \code{Psi = 0} uses the scale-invariant prior on \eqn{\Sigma}: \eqn{p(\Sigma) \propto |\Sigma|^{-(\nu+d+1)/2}}.
#'
#' The default value \code{nu = ncol(X)+1} for \code{kappa = 0} and \code{Psi = 0} makes \eqn{E[\mu|X]=\code{colMeans(X)}} and \eqn{E[\Sigma | X]=\code{var(X)}}.
#' @param nsamples the number of posterior draws.
#' @param X a data matrix with observations along the rows.
#' @param lambda location parameter. See Details.
#' @param kappa scale parameter. Defaults to \code{kappa = 0}.  See Details.
#' @param Psi scale matrix. Defaults to \code{Psi = 0}.  See Details.
#' @param nu degrees of freedom. Defaults to \code{nu = ncol(X)+1}.  See Details.
#' @return Returns a list with elements \code{mu} and \code{Sigma} of sizes \code{c(nsamples, length(lambda))} and \code{c(dim(Psi), nsamples)}.
#' @seealso \code{\link{rniw}}, \code{\link{niiw.post}}.
#' @example examples/niw.post.R
#' @export
niw.post <- function(nsamples, X, lambda, kappa, Psi, nu) {
  par <- niw.coeffs(X, lambda, kappa, Psi, nu)
  rniw(nsamples, par$lambda, par$kappa, par$Psi, par$nu)
}
