#' Point coordinates for a 2-D ellipse.
#'
#' Calculates coordinates of points for plotting a 2-dimensional ellipse based on user-defined parameters. Can be used for exploratory data analysis to produce ellipses at a given niche region size (e.g., \eqn{\alpha = 95\%}).
#'
#' @details This function provides the coordinates needed to plot a 2-dimensional ellipse based on user-defined parameters, such that \code{X = c(x,y)} satisfies the equation
#' \deqn{
#' (X-\mu)' V^{-1} (X-\mu) = C,
#' }
#' where \eqn{C=\code{qchisq(alpha, df = 2)}}.
#'
#' @param mu centre of ellipse. A vector of length 2.
#' @param V scale of ellipse. A 2x2 matrix. See Details.
#' @param alpha niche region size. See Details.
#' @param n number of points to return for plotting.
#' @return Returns a matrix of coordinates \code{cbind(x,y)} to plot a 2-dimensional ellipse.
#' @seealso \code{\link{niche.plot}}
#'
#' @example examples/ellipse.R
#' @export
ellipse <- function(mu, V, alpha = .95, n = 100) {
  tmp <- eigen(V)
  hlen <- sqrt(qchisq(alpha, df = 2)*tmp$val)
  theta <- atan2(tmp$vec[2,1], tmp$vec[1,1])
  t <- seq(0, 2*pi, len = n+1)
  x <- hlen[1] * cos(t)
  y <- hlen[2] * sin(t)
  alpha <- atan2(y, x)
  rad <- sqrt(x^2 + y^2)
  cbind(x = rad * cos(alpha + theta) + mu[1],
        y = rad * sin(alpha + theta) + mu[2])
}
