#' Plot for 2-d projection of niche regions.
#'
#' For one or more species, creates a series of plots: (i) the raw niche indicators (e.g., stable isotope) data, (ii) their density estimates, and (iii) 2-dimensional projections of probabilistic niche regions based on \eqn{n}-dimensionsional data.
#'
#' @details A set of plots is created for each pairwise combination of niche indicators.  Below the diagonal are scatterplots for each species, above the diagonal are ellipses corresponding to 2-d projections of the probabilistic niche regions.  The diagonal displays density estimates for each indicator, and optionally the raw 1-d data.  See Swanson et al. (2015) for detailed description of the probabilistic niche region.
#'
#' @param niche.par a list of length \code{nspecies}, each element of which in turn is a list with elements \code{mu} and \code{Sigma}.  Each of these will correspond to an ellipse being drawn for that species in the corresponding 2-d plane. See Example.
#'
#' @param niche.data a list of length \code{nspecies}, each element of which is a matrix with observations along the rows and niche indicators (e.g., stable isotopes) along the columns.
#' @param alpha size of the niche region to plot. Defaults to 0.95.
#' @param species.names names of the species. Defaults to \code{names(niche.par)}.
#' @param iso.names names of the niche indicators (or isotopes) to plot. Defaults to \code{colnames(niche.par)}.
#' @param lims two-row matrix of range limits for each niche indicator.  Defaults to include all ellipses.
#' @param col vector of colours in which each species will be drawn.
#' @param ndens number of points at which to evaluate kernel density estimates.
#' @param pfrac fraction of the plot on which to display 1-dimensional raw niche indicator data. \code{pfrac = 0} means don't display the raw data in 1-d.
#' @param xlab title of plot, located on the bottom.  Defaults to no title.
#' @return Returns a series of plots displaying niche indicator data and their probabilistic niche projections.
#' @references Swanson, H.K., Lysy, M., Stasko, A.D., Power, M., Johnson, J.D., and Reist, J.D. "A new probabilistic method for quantifying n-dimensional ecological niches and niche overlap." \emph{Ecology: Statistical Reports} 96:2 (2015): 318-324. \url{https://www.ncbi.nlm.nih.gov/pubmed/26240852}.
#' @seealso \code{\link{overlap.plot}}, \code{\link{niw.post}}, \code{\link{niiw.post}}.
#' @example examples/niche.plot.R
#' @export
niche.plot <- function(niche.par, niche.data, alpha = .95,
                       species.names, iso.names, lims,
                       col, ndens = 512, pfrac = 0, xlab) {
  niso <- ncol(niche.par[[1]]$mu)
  nspec <- length(niche.par)
  npts <- 100 # number of points for each ellipse
  nell <- sapply(niche.par, function(x) nrow(x$mu)) # number of ellipses per species
  if(missing(species.names)) species.names <- names(niche.par)
  if(missing(iso.names)) iso.names <- colnames(niche.par[[1]]$mu)
  # create all the ellipses to get the plot limits right.
  ell <- vector("list", nspec)
  names(ell) <- names(species.names)
  D <- combn(niso, 2)
  for(ii in 1:nspec) {
    ell.tmp <- array(NA, c(nell[ii], ncol(D), npts+1, 2))
    for(jj in 1:nell[ii]) {
      for(kk in 1:ncol(D)) {
        ell.coord <- ellipse(niche.par[[ii]]$mu[jj, D[,kk]],
                             V = niche.par[[ii]]$Sigma[D[,kk], D[,kk], jj],
                             alpha = alpha, n = npts)
        ell.tmp[jj,kk,,] <- ell.coord
      }
    }
    ell[[ii]] <- ell.tmp
  }
  # plot limits.
  if(missing(lims)) {
    # data
    lims <- array(sapply(niche.data, function(x) apply(x, 2, range)),
                  dim = c(2, niso, nspec))
    lims <- apply(lims, 2, range)
    # ellipse
    DD <- t(sapply(1:niso, function(ii) which(D == ii, arr.ind = TRUE)[1,]))
    elims <- array(sapply(ell, function(x) {
      tmp <- apply(x, c(4,2), range)
      rbind(as.matrix(tmp[1,,])[DD], as.matrix(tmp[2,,])[DD])
    }), dim = c(2, niso, nspec))
    elims <- apply(elims, 2, range)
    lims <- apply(rbind(lims, elims), 2, range)
  }
  # plots
  opar <- par(mfcol = c(niso,niso), mar = rep(.5, 4), oma = rep(4,4))
  for(ci in 1:niso) {
    for(ri in 1:niso) {
      # initialize plot
      plot.new()
      plot.window(lims[,ci], lims[,ri])
      if (ci == ri) {
        # diagonals: density plots
        xdens <- matrix(NA, ndens, nspec)
        ydens <- xdens
        for(ii in 1:nspec) {
          den <- density(niche.data[[ii]][,ci], n = ndens)
          xdens[,ii] <- den$x
          ydens[,ii] <- den$y
        }
        for(ii in 1:nspec) {
          ly <- par("usr")[1:2]
          ly[2] <- ly[1] + pfrac*(ly[2]-ly[1])
          ly[3] <- (ly[2]-ly[1])/nspec
          segments(x0 = niche.data[[ii]][,ci],
                   y0 = ly[1]+(ii-1)*ly[3], y1 = ly[1]+ii*ly[3], col = col[ii])
          ly <- ly[2] + ydens[,ii]/max(ydens)*(lims[2,ci]-ly[2])
          lines(xdens[,ii], ly, col = col[ii])
        }
      }
      if (ri > ci) {
        # lower triangle: point plots
        for(ii in 1:nspec) {
          points(niche.data[[ii]][,c(ci,ri)], col = col[ii], pch = 16)
        }
      }
      if (ri < ci) {
        # upper triangle: ellipses
        for(ii in 1:nspec) {
          for(jj in 1:nell[ii]) {
            lines(ell[[ii]][jj,which(D[1,] == ri & D[2,] == ci),,2:1], col = col[ii])
          }
        }
      }
      box()
      if(ci == niso) axis(side = 4) else axis(side = 4, labels = FALSE)
      if(ri == niso) axis(side = 1) else axis(side = 1, labels = FALSE)
      if(ri == 1) mtext(text = iso.names[ci], side = 3, line = 1)
      if(ci == 1) mtext(text = iso.names[ri], side = 2, line = 1)
    }
  }
  if(!missing(xlab)) {
    mtext(text = xlab, side = 1, outer = TRUE, line = 2.2, cex = .9)
  }
  legend(x = "topleft", legend = species.names, fill = col, bty = "n", cex = 1.25)
  par(opar) # reset par
}

# old version
## niche.plot <- function(niche.par, niche.data, alpha = .95,
##                        species.names, iso.names,
##                        col, ndens = 512, pfrac = 0, xlab) {
##   niso <- ncol(niche.par[[1]]$mu)
##   nspec <- length(niche.par)
##   npts <- 100 # number of points for each ellipse
##   nell <- sapply(niche.par, function(x) nrow(x$mu)) # number of ellipses per species
##   if(missing(species.names)) species.names <- names(niche.par)
##   if(missing(iso.names)) iso.names <- colnames(niche.par[[1]]$mu)
##   # create all the ellipses to get the plot limits right.
##   ell <- vector("list", nspec)
##   names(ell) <- names(species.names)
##   D <- combn(niso, 2)
##   for(ii in 1:nspec) {
##     ell.tmp <- array(NA, c(nell[ii], ncol(D), npts+1, 2))
##     for(jj in 1:nell[ii]) {
##       for(kk in 1:ncol(D)) {
##         ell.coord <- ellipse(niche.par[[ii]]$mu[jj, D[,kk]],
##                              V = niche.par[[ii]]$Sigma[D[,kk], D[,kk], jj],
##                              alpha = alpha, n = npts)
##         ell.tmp[jj,kk,,] <- ell.coord
##       }
##     }
##     ell[[ii]] <- ell.tmp
##   }
##   # plot limits.
##   lims <- array(sapply(niche.data, function(x) apply(x, 2, range)),
##                 dim = c(2, niso, nspec))
##   lims <- apply(lims, 2, range)
##   # plots
##   par(mfcol = c(niso,niso), mar = rep(.5, 4), oma = rep(4,4))
##   for(ci in 1:niso) {
##     for(ri in 1:niso) {
##       # initialize plot
##       plot.new()
##       plot.window(lims[,ci], lims[,ri])
##       if (ci == ri) {
##         # diagonals: density plots
##         xdens <- matrix(NA, ndens, nspec)
##         ydens <- xdens
##         for(ii in 1:nspec) {
##           den <- density(niche.data[[ii]][,ci], n = ndens)
##           xdens[,ii] <- den$x
##           ydens[,ii] <- den$y
##         }
##         for(ii in 1:nspec) {
##           ly <- par("usr")[1:2]
##           ly[2] <- ly[1] + pfrac*(ly[2]-ly[1])
##           ly[3] <- (ly[2]-ly[1])/nspec
##           segments(x0 = niche.data[[ii]][,ci],
##                    y0 = ly[1]+(ii-1)*ly[3], y1 = ly[1]+ii*ly[3], col = col[ii])
##           ly <- ly[2] + ydens[,ii]/max(ydens)*(lims[2,ci]-ly[2])
##           lines(xdens[,ii], ly, col = col[ii])
##         }
##       }
##       if (ri > ci) {
##         # lower triangle: point plots
##         for(ii in 1:nspec) {
##           points(niche.data[[ii]][,c(ci,ri)], col = col[ii], pch = 16)
##         }
##       }
##       if (ri < ci) {
##         # upper triangle: ellipses
##         for(ii in 1:nspec) {
##           for(jj in 1:nell[ii]) {
##             lines(ell[[ii]][jj,which(D[1,] == ri & D[2,] == ci),,2:1], col = col[ii])
##           }
##         }
##       }
##       box()
##       if(ci == niso) axis(side = 4) else axis(side = 4, labels = FALSE)
##       if(ri == niso) axis(side = 1) else axis(side = 1, labels = FALSE)
##       if(ri == 1) mtext(text = iso.names[ci], side = 3, line = 1)
##       if(ci == 1) mtext(text = iso.names[ri], side = 2, line = 1)
##     }
##   }
##   if(!missing(xlab)) {
##     mtext(text = xlab, side = 1, outer = TRUE, line = 2.2, cex = .9)
##   }
##   legend(x = "topleft", legend = species.names, fill = col, bty = "n", cex = 1.25)
## }
