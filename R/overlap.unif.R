#' Overlap calculation for uniform niche regions
#'
#' @details The overlap between niche regions \eqn{A} and \eqn{B} is defined as \eqn{vol(A \cap B)/vol(A \cup B)}, where the hypervolume of an n-dimensional region \eqn{S \in \mathbb R^n} is \eqn{vol(S) = \int_S dx}.  For elliptical niche regions, there are simple formulas for \eqn{vol(A)} and \eqn{vol(B)}.  Thus, we need only determine the volume of the intersection \eqn{vol(A \cap B)}, as the volume of the union is given by the formula \eqn{vol(A \cup B) = vol(A) + vol(B) - vol(A \cap B)}.
#'
#' For spherical niche regions, \eqn{vol(A \cap B)} has a closed-form expression (see References).  For elliptical regions, no such formula exists and a Monte Carlo method is used instead.  That is, \eqn{vol(A \cap B)} is calculated by sampling uniformly from \eqn{A}, then multiplying \eqn{vol(A)} by the fraction of sampled points which fall into \eqn{B}.
#'
#' While the uniform overlap metric is invariant to permutation of niche regions \eqn{A} and \eqn{B}, the accuracy of the Monte Carlo calculation of \eqn{vol(A \cap B)} is not: higher accuracy is obtained when a higher fraction of sampled points are in the opposite niche region.  \code{overlap.unif} does not attempt to determine for which region this is the case, though the choice can be informed by plotting the niche regions, e.g., with \code{\link{niche.plot}}.
#' @param muA,muB mean of niche regions.
#' @param SigmaA,SigmaB variance matrix of elliptical niche regions.
#' @param alphaA,alphaB probabilistic size of niche regions.
#' @param nprob number of uniform draws from niche region \code{A}.
#' @return A Monte Carlo estimate of the niche overlap.
#' @name overlap.unif
#' @examples
#' # spherical case: compare Monte Carlo method to analytic formula
#'
#' d <- 2 # 2D example
#' mA <- rnorm(d)
#' mB <- rnorm(d)
#' sigA <- rexp(1)
#' SigA <- sigA^2 * diag(d)
#' sigB <- rexp(1)
#' SigB <- sigB^2 * diag(d)
#'
#' # plot circles
#' ellA <- ellipse(mA, SigA)
#' ellB <- ellipse(mB, SigB)
#' plot(0, type = "n",
#'      xlim = range(ellA[,1], ellB[,1]),
#'      ylim = range(ellA[,2], ellB[,2]), xlab = "x", ylab = "y")
#' lines(ellA, col = "red")
#' lines(ellB, col = "blue")
#' legend("topright", legend = c("niche A", "niche B"),
#'        fill = c("red", "blue"), bg = "white")
#'
#' # compare niche calculations
#' overlap.sphere(mA, sigA, mB, sigB)
#' overlap.unif(mA, SigA, mB, SigB, nprob = 1e5)
#' @export
overlap.unif <- function(muA, SigmaA, muB, SigmaB, alphaA = .95, alphaB = .95,
                         nprob) {
  # niche sizes
  volA <- niche.size(Sigma = SigmaA, alpha = alphaA)
  volB <- niche.size(Sigma = SigmaB, alpha = alphaB)
  # intersection: simulate data from one niche
  simA <- niche.runif(nprob, mu = muA, Sigma = SigmaA, alpha = alphaA)
  # fraction of pts in the other
  Z <- backsolve(chol(SigmaB), t(simA) - muB, transpose = TRUE)
  propAB <- mean(colSums(Z*Z) < qchisq(p = alphaB, df = length(muB)))
  # overlap calculation
  volAB <- propAB * volA
  volAB / (volA + volB - volAB)
}

#' @rdname overlap.unif
#' @param sigmaA,sigmaB standard deviations (scalars) of spherical niche regions.
#' @references Li, S. "Concise formulas for the area and volume of a hyperspherical cap." \emph{Asian Journal of Mathematics & Statistics} 4.1 (2011): 66-70. \url{http://dx.doi.org/10.3923/ajms.2011.66.70}.
#' @export
overlap.sphere <- function(muA, sigmaA, muB, sigmaB,
                           alphaA = .95, alphaB = .95) {
  n <- length(muA) # number of dimensions
  # niche sizes
  volA <- niche.size(Sigma = sigmaA^2 * diag(n), alpha = alphaA)
  volB <- niche.size(Sigma = sigmaB^2 * diag(n), alpha = alphaB)
  # volume of intersection
  rA <- sigmaA * sqrt(qchisq(alphaA, df = n)) # sphere radii
  rB <- sigmaB * sqrt(qchisq(alphaB, df = n))
  dd <- sqrt(sum((muA - muB)^2)) # distance between sphere centers
  if(dd > rA + rB) {
    volAB <- 0
  } else if(dd < abs(rA-rB)) {
    volAB <- min(volA, volB)
  } else {
    cA <- .5 * (dd^2 + rA^2 - rB^2)/dd
    cB <- .5 * (dd^2 - rA^2 + rB^2)/dd
    volAB <- vol.cap(rA, cA, n) + vol.cap(rB, cB, n)
  }
  volAB / (volA + volB - volAB)
}

# volume of an n-dimensional sphere cap
vol.cap <- function(r, a, n) {
  vol <- .5 * pbeta((a/r)^2, .5, .5*(n+1), lower.tail = FALSE)
  if(a < 0) vol <- 1 - vol
  (sqrt(pi) * r)^n/gamma(1+.5*n) * vol
}

# cseg <- function(R, d) R^2 * acos(d/R) - d*sqrt(R^2 - d^2)

## # 1D example
## d <- 1
## mA <- rnorm(d)
## mB <- rnorm(d)
## sigA <- rexp(1)
## SigA <- sigA^2 * diag(d)
## sigB <- rexp(1)
## SigB <- sigB^2 * diag(d)
## nicheA <- mA + c(-1,1) * sigA * sqrt(qchisq(.95, 1))
## nicheB <- mB + c(-1,1) * sigB * sqrt(qchisq(.95, 1))

## plot(0, type = "n", xlim = range(nicheA, nicheB), ylim = c(0,3),
##      xlab = "x", ylab = "y")
## lines(nicheA, c(1,1), col = "red", lwd = 3)
## lines(nicheB, c(2,2), col = "blue", lwd = 3)
## points(c(mA, mB), c(1,2), pch = 16, cex = 2)

## overlap.sphere(mA, sigA, mB, sigB)
## overlap.unif(mA, SigA, mB, SigB, nprob = 1e5)

## # 3d example

## require(rgl)
## open3d()

## d <- 3
## mA <- rnorm(d)
## ## mB <- mA + .1
## mB <- rnorm(d)
## sigA <- rexp(1)
## SigA <- sigA^2 * diag(d)
## ## sigB <- sigA
## sigB <- rexp(1)
## SigB <- sigB^2 * diag(d)
## clear3d()
## # open3d()
## shade3d(ellipse3d(x = SigA, centre = mA), col = "yellow")
## shade3d(ellipse3d(x = SigB, centre = mB), col = "green")

## overlap.sphere(mA, sigA, mB, sigB)
## overlap.unif(mA, SigA, mB, SigB, nprob = 1e5)
