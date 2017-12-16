\name{niche.par.plot}
\alias{niche.par.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
niche.par.plot(niche.par, plot.mu = TRUE, plot.Sigma = TRUE, plot.index, col, ndens = 512, ylab)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{niche.par}{
%%     ~~Describe \code{niche.par} here~~
}
  \item{plot.mu}{
%%     ~~Describe \code{plot.mu} here~~
}
  \item{plot.Sigma}{
%%     ~~Describe \code{plot.Sigma} here~~
}
  \item{plot.index}{
%%     ~~Describe \code{plot.index} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{ndens}{
%%     ~~Describe \code{ndens} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (niche.par, plot.mu = TRUE, plot.Sigma = TRUE, plot.index, 
    col, ndens = 512, ylab) 
{
    niso <- ncol(niche.par[[1]]$mu)
    nsmp <- length(niche.par)
    if (!missing(plot.index)) {
        if (length(plot.index) == 1) {
            plot.mu <- TRUE
            plot.Sigma <- FALSE
        }
        else if (length(plot.index) == 2) {
            plot.mu <- FALSE
            plot.Sigma <- TRUE
        }
        else {
            stop("Incorrect specification of plot.index.  Must be a numeric vector of length 1 (plot.mu) or 2 (plot.Sigma).")
        }
        nc <- 1
        nr <- 1
    }
    else {
        nc <- niso
        nr <- 0
        if (plot.mu) 
            nr <- nr + 1
        if (plot.Sigma) 
            nr <- nr + niso
    }
    if (missing(ylab)) {
        if (plot.mu) {
            ylab.mu <- sapply(1:niso, function(ii) parse(text = paste0("p(mu[", 
                ii, "]*\" | \"*X)")))
        }
        if (plot.Sigma) {
            ylab.Sigma <- apply(as.matrix(expand.grid(1:niso, 
                1:niso))[, 2:1], 1, function(II) {
                paste0(II[1], "*", II[2])
            })
            ylab.Sigma <- sapply(ylab.Sigma, function(ii) {
                parse(text = paste0("p(Sigma[", ii, "]*\" | \"*X)"))
            })
        }
    }
    else {
        ylab <- rep(ylab, len = niso * (niso + 1))
        ylab.mu <- ylab[1:niso]
        ylab.Sigma <- ylab[-(1:niso)]
    }
    densx <- matrix(NA, ndens, nsmp)
    densy <- densx
    if ((nr > 1) || (nc > 1)) 
        par(mfrow = c(nr, nc))
    if (plot.mu) {
        for (ii in 1:niso) {
            if (missing(plot.index) || plot.index == ii) {
                for (kk in 1:nsmp) {
                  dens <- density(niche.par[[kk]]$mu[, ii], n = ndens)
                  densx[, kk] <- dens$x
                  densy[, kk] <- dens$y
                }
                plot(densx, densy, type = "n", xlab = parse(text = paste0("mu[", 
                  ii, "]")), ylab = ylab.mu[ii])
                for (kk in 1:nsmp) lines(densx[, kk], densy[, 
                  kk], col = col[kk])
            }
        }
    }
    if (plot.Sigma) {
        for (ii in 1:niso) {
            for (jj in 1:niso) {
                if (missing(plot.index) || all(plot.index == 
                  c(ii, jj))) {
                  for (kk in 1:nsmp) {
                    dens <- density(niche.par[[kk]]$Sigma[ii, 
                      jj, ], n = ndens)
                    densx[, kk] <- dens$x
                    densy[, kk] <- dens$y
                  }
                  plot(densx, densy, type = "n", xlab = parse(text = paste0("Sigma[", 
                    ii, "*", jj, "]")), ylab = ylab.Sigma[(ii - 
                    1) * niso + jj])
                  for (kk in 1:nsmp) lines(densx[, kk], densy[, 
                    kk], col = col[kk])
                }
            }
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
