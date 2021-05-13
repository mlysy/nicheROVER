#' Random draws from the posterior distribution with Normal-Independent-Inverse-Wishart (NIIW) prior.
#'
#' Given iid \eqn{d}-dimensional niche indicators  \eqn{X = (X_1,\ldots,X_N)} with \eqn{X_i \sim N(\mu, \Sigma)}, this function generates random draws from \eqn{p(\mu,\Sigma | X)} for the Normal-Independent-Inverse-Wishart (NIIW) prior.
#'
#' @details The NIIW distribution \eqn{p(\mu, \Sigma | \lambda, \kappa, \Psi, \nu)} is defined as
#' \deqn{
#' \Sigma \sim W^{-1}(\Psi, \nu), \quad \mu | \Sigma \sim N(\lambda, \Omega).
#' }
#' The default value `Omega = 0` uses the Lebesque prior on \eqn{\mu}: \eqn{p(\mu) \propto 1}.  In this case the NIW and NIIW priors produce identical resuls, but [niw.post()] is faster.
#'
#' The default value `Psi = 0` uses the scale-invariant prior on \eqn{\Sigma}: \eqn{p(\Sigma) \propto |\Sigma|^{-(\nu+d+1)/2}}.
#'
#' The default value `nu = ncol(X)+1` for `Omega = 0` and `Psi = 0` makes \eqn{E[\mu|X]=`colMeans(X)`} and \eqn{E[\Sigma | X]=`var(X)`}.
#'
#' Random draws are obtained by a Markov chain Monte Carlo (MCMC) algorithm; specifically, a Gibbs sampler alternates between draws from \eqn{p(\mu | \Sigma, X)} and \eqn{p(\Sigma | \mu, X)}, which are Normal and Inverse-Wishart distributions respectively.
#'
#' @param nsamples The number of posterior draws.
#' @param X A data matrix with observations along the rows.
#' @param lambda Mean of \eqn{\mu}. See 'Details'.
#' @param Omega Variance of \eqn{\mu}. Defaults to `Omega = 0`.  See 'Details'.
#' @param Psi Scale matrix of \eqn{\Sigma}. Defaults to `Psi = 0`.  See 'Details'.
#' @param nu Degrees of freedom of \eqn{\Sigma}. Defaults to `nu = ncol(X)+1`.  See 'Details'.
#' @param mu0 Initial value of \eqn{\mu} to start the Gibbs sampler.  See 'Details'.
#' @param burn Burn-in for the MCMC sampling algorithm.  Either an integer giving the number of initial samples to discard, or a fraction with `0 < burn < 1`.  Defaults to `burn = floor(nsamples/10)`.
#' @return Returns a list with elements `mu` and `Sigma` of sizes `c(nsamples, length(lambda))` and `c(dim(Psi), nsamples)`.
#' @seealso [niw.post()], [rwish()].
#' @example examples/niiw.post.R
#' @export
niiw.post <- function(nsamples, X, lambda, Omega, Psi, nu, mu0 = lambda, burn) {
  # sufficient statistics
  X <- as.matrix(X)
  d <- ncol(X)
  N <- nrow(X)
  Xbar <- colMeans(X)
  S <- t(X)-Xbar
  S <- S %*% t(S)
  # local variables
  mu <- rep(0, d)
  Sigma <- matrix(0, d, d)
  if(missing(burn)) burn <- .1
  if(burn < 1) burn <- floor(nsamples*burn)
  mu <- mu0
  # output
  mu.out <- matrix(NA, nsamples, d)
  Sigma.out <- array(NA, dim = c(d, d, nsamples))
  # main for loop
  for(ii in (-burn+1):nsamples) {
    # sample Sigma
    Psi2 <- N * ((mu-Xbar) %o% (mu-Xbar)) + S + Psi
    nu2 <- N+nu
    Sigma <- matrix(rwish(1, Psi2, nu2, inv = TRUE), d, d)
    # sample mu
    Sigma2 <- Sigma/N
    if(!all(Omega == 0)) {
      B <- Sigma2 %*% solve(Omega + Sigma2)
    } else B <- matrix(0, d, d)
    IB <- diag(d)-B
    lambda2 <- c(IB %*% Xbar)
    if(!all(Omega == 0)) lambda2 <- lambda2 + c(B %*% lambda)
    Omega2 <- IB %*% Sigma2
    mu <- .rmvn(mean = lambda2, sigma = Omega2)
    # store
    if(ii > 0) {
      mu.out[ii,] <- mu
      Sigma.out[,,ii] <- Sigma
    }
  }
  list(mu = mu.out, Sigma = Sigma.out)
}
