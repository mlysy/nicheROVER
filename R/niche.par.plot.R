#' Plot for niche parameters.
#'
#' For one or more species, plots some or all of the niche parameters \eqn{\mu} and \eqn{\Sigma}.
#'
#' @param niche.par list with \code{nspecies = length(niche.par)}, each element of which is a list with parameters \code{mu} and \code{Sigma}.  See Details.
#' @param plot.mu logical.  If \code{TRUE}, plot the distribution of \eqn{\mu} for each niche indicator (e.g., stable isotope).  See Details.
#' @param plot.Sigma logical.  If \code{TRUE}, plot the distribution of \eqn{\Sigma} for each niche indicator.  See Details.
#' @param plot.index either a scalar of a numeric vector of length 2.  If \code{plot.index = i} then plot the distribution of \eqn{\mu_i}.  If \code{plot.index = c(i,j)} then plot the distribution of \eqn{\Sigma_{ij}}.
#' @param col vector of colors in which to plot each species.
#' @param ndens number of points at which to evaluate density estimates.
#' @param ylab optional label for \eqn{y}-axis.  If missing, defaults to \eqn{p(\mu_i | X)} and \eqn{p(\Sigma_{ij} | X)}.
#'@return Returns a plot of the distribution of some or all niche parameters.
#'
#' @details \code{niche.par} is a list, each element of which is a distribution of niche parameters.  That is, \code{names(niche.par[[1]]) = c("mu", "Sigma")}, and if \code{niso} is the number of niche indicators (e.g., stable isotopes), then \code{dim(niche.par[[1]]$mu) = c(nsamples, niso)} and \code{dim(niche.par[[1]]$Sigma) = c(niso, niso, nsamples)}.
#'
#' @seealso \code{\link{niw.post}}, \code{\link{niiw.post}} for niche parameter output, \code{density} in the \code{R} \code{base} package for density estimation from sample data.
#'
#' @example examples/niche.par.plot.R
#' @export
niche.par.plot <- function(niche.par, plot.mu = TRUE, plot.Sigma = TRUE, plot.index,
                           col, ndens = 512, ylab) {
  niso <- ncol(niche.par[[1]]$mu)
  nsmp <- length(niche.par)
  # determine the number of rows and columns in plot
  if(!missing(plot.index)) {
    if(length(plot.index) == 1) {
      plot.mu <- TRUE
      plot.Sigma <- FALSE
    } else if(length(plot.index) == 2) {
      plot.mu <- FALSE
      plot.Sigma <- TRUE
    } else {
      stop("Incorrect specification of plot.index.  Must be a numeric vector of length 1 (plot.mu) or 2 (plot.Sigma).")
    }
    nc <- 1
    nr <- 1
  } else {
    nc <- niso
    nr <- 0
    if(plot.mu) nr <- nr + 1
    if(plot.Sigma) nr <- nr + niso
  }
  # determine ylab
  if(missing(ylab)) {
    if(plot.mu) {
      ylab.mu <- sapply(1:niso, function(ii)
                        parse(text = paste0("p(mu[", ii, "]*\" | \"*X)")))
    }
    if(plot.Sigma) {
      ylab.Sigma <- apply(as.matrix(expand.grid(1:niso, 1:niso))[,2:1], 1, function(II) {
        paste0(II[1], "*", II[2])
      })
      ylab.Sigma <- sapply(ylab.Sigma, function(ii) {
        parse(text = paste0("p(Sigma[", ii, "]*\" | \"*X)"))
      })
    }
  } else {
    ylab <- rep(ylab, len = niso*(niso+1))
    ylab.mu <- ylab[1:niso]
    ylab.Sigma <- ylab[-(1:niso)]
  }
  # plot
  densx <- matrix(NA, ndens, nsmp)
  densy <- densx
  if((nr > 1) || (nc > 1)) par(mfrow = c(nr, nc))
  # mu
  if(plot.mu) {
    for(ii in 1:niso) {
      if(missing(plot.index) || plot.index == ii) {
        for(kk in 1:nsmp) {
          dens <- density(niche.par[[kk]]$mu[,ii], n = ndens)
          densx[,kk] <- dens$x
          densy[,kk] <- dens$y
        }
        plot(densx, densy, type = "n",
             xlab = parse(text = paste0("mu[", ii, "]")),
             ylab = ylab.mu[ii])
        for(kk in 1:nsmp) lines(densx[,kk], densy[,kk], col = col[kk])
      }
    }
  }
  if(plot.Sigma) {
    for(ii in 1:niso) {
      for(jj in 1:niso) {
        if(missing(plot.index) || all(plot.index == c(ii,jj))) {
          for(kk in 1:nsmp) {
            dens <- density(niche.par[[kk]]$Sigma[ii,jj,], n = ndens)
            densx[,kk] <- dens$x
            densy[,kk] <- dens$y
          }
          plot(densx, densy, type = "n",
               xlab = parse(text = paste0("Sigma[", ii, "*", jj, "]")),
               ylab = ylab.Sigma[(ii-1)*niso + jj])
          for(kk in 1:nsmp) lines(densx[,kk], densy[,kk], col = col[kk])
        }
      }
    }
  }
}
