#'@title Random draws from the posterior distribution with Normal-Independent-Inverse-Wishart (NIIW) prior.
#'
#'@description Given iid \eqn{d}-dimensional niche indicators  \eqn{X = (X_1,\ldots,X_N)} with \eqn{X_i \sim N(\mu, \Sigma)},
#'this function generates random draws from \eqn{p(\mu,\Sigma | X)} for the Normal-Independent-Inverse-Wishart (NIIW) prior.
#'@details The NIIW distribution \eqn{p(\mu, \Sigma | \lambda, \kappa, \Psi, \nu)} is defined as
#'\deqn{\Sigma \sim W^{-1}(\Psi, \nu), \quad \mu | \Sigma \sim N(\lambda, \Omega).}
#'The default value \code{Omega = 0} uses the Lebesque prior on \eqn{\mu}: \eqn{p(\mu) \propto 1}.  In this case the NIW and NIIW priors produce identical resuls, but \code{\link{niw.post}} is faster.
#'The default value \code{Psi = 0} uses the scale-invariant prior on \eqn{\Sigma}: \eqn{p(\Sigma) \propto |\Sigma|^{-(\nu+d+1)/2}}.
#'The default value \code{nu = ncol(X)+1} for \code{Omega = 0} and \code{Psi = 0} makes \eqn{E[\mu|X]=\code{colMeans(X)}} and \eqn{E[\Sigma | X]=\code{var(X)}}.
#'Random draws are obtained by a Markov chain Monte Carlo (MCMC) algorithm; specifically,
#'a Gibbs sampler alternates between draws from \eqn{p(\mu | \Sigma, X)} and \eqn{p(\Sigma | \mu, X)}, which are Normal and Inverse-Wishart distributions respectively.
#'@param nsamples the number of posterior draws.
#'@param X a data matrix with observations along the rows.
#'@param lambda mean of mu. See Details.
#'@param Omega variance of mu. Defaults to \code{Omega = 0}.  See Details.
#'@param Psi scale matrix of Sigma. Defaults to \code{Psi = 0}.  See Details.
#'@param nu degrees of freedom of Sigma. Defaults to \code{nu = ncol(X)+1}.  See Details.
#'@param mu0 initial value of mu to start the Gibbs sampler.  See Details.
#'@param burn burn-in for the MCMC sampling algorithm.  Either an integer giving the number of initial samples to discard, or a fraction with \code{0 < burn < 1}.  Defaults to \code{burn = floor(nsamples/10)}.
#'@return Returns a list with elements \code{mu} and \code{Sigma} of sizes \code{c(nsamples, length(lambda))} and \code{c(dim(Psi), nsamples)}.
#'@seealso \code{\link{niw.post}}, \code{\link{rwish}}.
#'@examples
#'# simulate normal data with mean and variance (mu0, Sigma0)
#'d <- 4
#'mu0 <- rnorm(d)
#'Sigma0 <- matrix(rnorm(d^2), d, d)
#'Sigma0 <- Sigma0 %*% t(Sigma0)
#'N <- 1e2
#'X <- matrix(rnorm(N*d), N, d) # iid N(0,1)
#'X <- t(t(X %*% chol(Sigma0)) + mu0) # each row is N(mu0, Sigma)
#'
#'# prior parameters
#'# flat prior on mu
#'lambda <- 0
#'Omega <- 0
#'# informative prior on Sigma
#'Psi <- crossprod(matrix(rnorm(d^2), d, d))
#'nu <- 5
#'
#'# sample from NIIW posterior
#'nsamples <- 2e3
#'system.time({
#'  siiw <- niiw.post(nsamples, X, lambda, Omega, Psi, nu, burn = 100)
#'})
#'
#'# sample from NIW posterior
#'kappa <- 0
#'system.time({
#'  siw <- niw.post(nsamples, X, lambda, kappa, Psi, nu)
#'})
#'
#'# check that posteriors are the same
#'
#'# p(mu | X)
#'clrs <- c("black", "red")
#'par(mar = c(4.2, 4.2, 2, 1)+.1)
#'niche.par.plot(list(siiw, siw), col = clrs, plot.mu = TRUE, plot.Sigma = FALSE)
#'legend(x = "topright", legend = c("NIIW Prior", "NIW Prior"), fill = clrs)
#'
#'# p(Sigma | X)
#'par(mar = c(4.2, 4.2, 2, 1)+.1)
#'niche.par.plot(list(siiw, siw), col = clrs, plot.mu = FALSE, plot.Sigma = TRUE)
#'legend(x = "topright", legend = c("NIIW Prior", "NIW Prior"), fill = clrs)
#'@export
niiw.post <- function(nsamples, X, lambda, Omega, Psi, nu, mu0 = lambda, burn) {
  # sufficient statistics
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
