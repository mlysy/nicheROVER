#' Calculate the size of an elliptical niche region.
#'
#' @details For a given niche region \eqn{N_R}, the niche size is defined as the hypervolume of this region: \eqn{N_S = \int_{x \in N_R} d x}.
#' @param Sigma variance matrix for normally distributed niche axes.
#' @param alpha probabilistic niche size.
#' @examples
#' # for each species, size of 95% niche region using sample variance
#' tapply(1:nrow(fish), fish$species, function(ind) {
#'   X <- fish[ind,2:4] # all measurements for given species
#'   Sighat <- var(X) # sample variance
#'   niche.size(Sigma = Sighat, alpha = .95)
#' })
#' @return hypervolume niche size.
#' @export
niche.size <- function(Sigma, alpha = .95) {
  n <- nrow(Sigma)
  sz <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  sz <- .5 * (sz + n * (log(pi) + log(qchisq(alpha, df = n)))) - lgamma(.5*n+1)
  exp(sz)
}
