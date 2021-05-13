#' Random draws from a Wishart (or Inverse-Wishart) distribution.
#'
#' Generates a random samples from a Wishart distribution defined as \eqn{W(\Psi, \nu)}, or an Inverse-Wishart distribution defined as \eqn{W^{-1}(\Psi, \nu)}.
#'
#' @details Setting `inv = TRUE` replaces \eqn{\Psi} by \eqn{Psi^{-1}} and inverts the output random matrices, such that they are being generated from an Inverse-Wishart \eqn{W^{-1}(\Psi, \nu)} distribution.
#'
#' @param n number of samples to draw.
#' @param Psi scale matrix.
#' @param nu degrees of freedom.
#' @param inv logical. Setting `inv = TRUE` returns random matrices from an Inverse-Wishart distribution. See Details.
#' @seealso [rniw()]
#' @return Returns an array of Wishart (or Inverse-Wishart) draws of size `c(nrow(Psi),ncol(Psi),n)`.
#' @example examples/rwish.R
#' @export
rwish <- function(n, Psi, nu, inv = FALSE) {
  if(inv) Psi <- solve(Psi)
  U <- chol(Psi)
  d <- nrow(Psi)
  ans <- array(0, dim = c(d, d, n))
  if(!is.null(dimnames(Psi))) dimnames(ans) <- c(dimnames(Psi), list(NULL))
  ans[rep(upper.tri(Psi), n)] <- rnorm(n*d*(d-1)/2)
  ans[rep(!lower.tri(Psi, diag = FALSE) &
            !upper.tri(Psi, diag = FALSE), n)] <- sqrt(rchisq(n*d, df = nu-1:d+1))
  for(ii in 1:n) {
    tmp <- ans[,,ii] %*% U
    if(inv) tmp <- backsolve(tmp, diag(d), transpose = TRUE)
    ans[,,ii] <- crossprod(tmp)
  }
  ans
}
