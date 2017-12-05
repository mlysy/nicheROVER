#--- todo for nicheROVER version 2.0 --------------------------------------------

# mlysy@uwaterloo.ca, july 2016

# calculate niche size N_S, defined as \eq{\int_{x \in N_R} d x}
# Sigma: variance matrix for normally distributed niche axes
# alpha: probabilistic niche size
niche.size <- function(Sigma, alpha = .95) {
  n <- nrow(Sigma)
  sz <- as.numeric(determinant(Sigma, log = TRUE)$modulus)
  sz <- .5 * (sz + n * (log(pi) + log(qchisq(alpha, df = n)))) - lgamma(.5*n+1)
  exp(sz)
}

#--- improved niche.plot --------------------------------------------------------

# 1. resets the "par" after use
# 2. custom limits, or use ellipses

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
