#' Posterior coefficients of the Normal-Inverse-Wishart distribution with its conjugate prior.
#'
#' Given iid \eqn{d}-dimensional niche indicators \eqn{X = (X_1,\ldots,X_N)} with \eqn{X_i \sim N(\mu, \Sigma)}, this function calculates the coefficients of the Normal-Inverse-Wishart (NIW) posterior \eqn{p(\mu, \Sigma | X)} for a conjugate NIW prior.  Together with [niw.mom()], this can be used to rapidly compute the point estimates \eqn{E[\mu | X]} and \eqn{E[\Sigma | X]}.
#'
#' @details The NIW distribution \eqn{p(\mu, \Sigma | \lambda, \kappa, \Psi, \nu)} is defined as
#' \deqn{
#' \Sigma \sim W^{-1}(\Psi, \nu), \quad \mu | \Sigma \sim N(\lambda, \Sigma/\kappa).
#' }
#' The default value `kappa = 0` uses the Lebesque prior on \eqn{\mu}: \eqn{p(\mu) \propto 1}.
#'
#' The default value `Psi = 0` uses the scale-invariant prior on \eqn{\Sigma}: \eqn{p(\Sigma) \propto |\Sigma|^{-(\nu+d+1)/2}}.
#'
#' The default value `nu = ncol(X)+1` for `kappa = 0` and `Psi = 0` makes \eqn{E[\mu|X]=`colMeans(X)`} and \eqn{E[\Sigma | X]=`var(X)`}.
#'
#' @param X a data matrix with observations along the rows.
#' @param lambda location parameter. See Details.
#' @param kappa scale parameter. Defaults to `kappa = 0`.  See Details.
#' @param Psi scale matrix. Defaults to `Psi = 0`.  See Details.
#' @param nu degrees of freedom. Defaults to `nu = ncol(X)+1`.  See Details.
#'
#' @return Returns a list with elements `lambda`, `kappa`, `Psi`, `nu` corresponding to the coefficients of the NIW posterior distribution \eqn{p(\mu, \Sigma | X)}.
#' @seealso [rniw()], [niw.mom()], [niw.post()].
#' @example examples/niw.coeffs.R
#' @export
niw.coeffs <- function(X, lambda, kappa, Psi, nu) {
  X <- as.matrix(X)
  d <- ncol(X)
  N <- nrow(X)
  if(missing(kappa)) kappa <- 0
  if(missing(Psi)) {
    if(N <= d) {
      stop("Must have more observations than niche dimensions when prior parameter Psi is missing.")
    }
    Psi <- 0
  }
  if(missing(nu)) nu <- ncol(X)+1
  # sufficient statistics
  Xbar <- colMeans(X)
  S <- t(X)-Xbar
  S <- S %*% t(S)
  # posterior parameter values
  Psi2 <- Psi + S
  lambda2 <- N*Xbar
  if(kappa != 0) {
    Psi2 <- Psi2 + (N*kappa)/(N+kappa) * (Xbar-lambda) %*% t(Xbar-lambda)
    lambda2 <- lambda2 + kappa*lambda
  }
  lambda2 <- lambda2/(N+kappa)
  nu2 <- N+nu-(kappa==0)
  kappa2 <- N+kappa
  list(lambda = lambda2, kappa = kappa2, Psi = Psi2, nu = nu2)
}
