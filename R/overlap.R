#' Monte Carlo calculation of niche region overlap metrics.
#'
#' Calculates the distribution of a niche region overlap metric for each pairwise species combination and user-specified niche region sizes.
#'
#' @details The overlap metric is the probability that a randomly drawn individual from species \eqn{A} will be found within the niche region of species \eqn{B} (for a given niche region size, e.g., \code{alpha = .95}).  It is a single number which is a function of the parameters for each species, \eqn{\Theta_A = (\mu_A, \Sigma_A)} and \eqn{\Theta_B = (\mu_B, \Sigma_B)}.  This number is difficult to calculate directly, but easy to approximate stochastically by generating \code{nprob} draws from the distribution of species \eqn{A} and counting the fraction of them which fall in the niche region of species \eqn{B}.
#'
#' Typically the true values of \eqn{\Theta_A} and \eqn{\Theta_B} are unknown and must be estimated from the data. Thus, the overlap metric is calculated for \code{nreps} combinations of samples from \eqn{p(\Theta_A | X)} and \eqn{p(\Theta_B | X)} which are supplied in \code{niche.par}.
#'
#' See Swanson et al. (2015) for a detailed description of niche overlap and its calculation.
#'
#' @param niche.par a list with \code{nspecies = length(niche.par)}, each element of which in turn is a list with elements \code{mu} and \code{Sigma}.  See Details.
#' @param nreps the number of overlap metric calculations for each species.  Defaults to the smallest number of parameter samples supplied by \code{niche.par}.  See Details.
#' @param nprob the number of normal draws for each Monte Carlo overlap metric calculation.  See Details.
#' @param alpha scalar or vector of niche region sizes for calculating the niche overlap metric. Defaults to 0.95.
#' @param species.names names of the species. Defaults to \code{names(niche.par)}.
#' @param norm.redraw logical. If \code{FALSE}, the same \code{nprob*nspecies} iid \eqn{N(0,1)} draws are used for each calculation of the overlap metric. This increases the Monte Carlo error, but the procedure is about 1.5x faster.  Defaults to \code{TRUE}.
#' @return Returns an array of size \code{c(nspecies, nspecies, nreps, nlevels)}, where \code{nlevels} is the number of alpha levels at which to calculate the overlap metric.  For each of the last two dimensions of the output array, the first two dimensions form an \code{nspecies} by \code{nspecies} matrix giving each pairwise calculation of overlap metric between two species for given \eqn{\Theta_A}, \eqn{\Theta_B}, and \code{alpha}. In each of these matrices, Species \eqn{A} is along the rows of this matrix and Species \eqn{B} is along the columns.
#'
#' @references Swanson, H.K., Lysy, M., Stasko, A.D., Power, M., Johnson, J.D., and Reist, J.D. "A new probabilistic method for quantifying n-dimensional ecological niches and niche overlap." \emph{Ecology: Statistical Reports} 96:2 (2015): 318-324. \url{https://www.ncbi.nlm.nih.gov/pubmed/26240852}.
#' @seealso \code{\link{overlap.plot}}, \code{\link{niw.post}}, \code{\link{niiw.post}}.
#' @example examples/overlap.R
#' @export
overlap <- function(niche.par, nreps, nprob, alpha = 0.95,
                    species.names, norm.redraw = TRUE) {
  niso <- ncol(niche.par[[1]]$mu) # number of isotopes
  nspec <- length(niche.par) # number of species
  nlevels <- length(alpha) # number of levels
  if(missing(species.names)) species.names <- names(niche.par)
  # temporary variables
  mu <- matrix(NA, niso, nspec)
  x <- array(NA, dim = c(niso, nprob, nspec))
  C <- array(NA, dim = c(niso, niso, nspec))
  Zsq <- array(NA, dim = c(nprob, nspec-1, nlevels))
  # constants
  qlevels <-  qchisq(alpha, df = niso)
  qlevels <- array(rep(qlevels, each = nprob*(nspec-1)),
                   dim = c(nprob, nspec-1, nlevels))
  if(!norm.redraw) z <- matrix(rnorm(niso*nprob), niso, nprob)
  # output
  over <- array(NA, dim = c(nspec, nspec, nreps, nlevels))
  dimnames(over) <- list(species.names, species.names, NULL,
                         paste0(alpha*100, "%"))
  names(dimnames(over)) <- c("Species A", "Species B", "", "alpha")
  # subsample the parameters posteriors of each niche nreps times
  ind <- sapply(niche.par, function(mv) {
    sample(nrow(mv$mu), nreps, replace = TRUE)
  })
  # calculate each overlap
  for(jj in 1:nreps) {
    # get means, t(chol(variance)), and draws for each species
    for(ii in 1:nspec) {
      mu[,ii] <- niche.par[[ii]]$mu[ind[jj,ii],]
      C[,,ii] <- t(chol(niche.par[[ii]]$Sigma[,,ind[jj,ii]]))
      if(norm.redraw) z <- matrix(rnorm(niso*nprob), niso, nprob)
      x[,,ii] <- C[,,ii] %*% z + mu[,ii]
    }
    # check whether the draw of each species is within the ellipse of every other
    for(ii in 1:nspec) {
      Z <- backsolve(C[,,ii], matrix(x[,,-ii], nrow = niso)-mu[,ii], upper.tri = FALSE)
      Zsq[] <- colSums(Z^2)
      over[-ii,ii,jj,] <- colMeans(Zsq < qlevels)
    }
  }
  if(dim(over)[4] == 1) over <- over[,,,1]
  over
}
