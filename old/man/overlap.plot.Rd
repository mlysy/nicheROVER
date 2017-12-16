\name{overlap.plot}
\alias{overlap.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
overlap.plot(over.stat, nbreaks = 50, equal.axis = FALSE, species.names, col, mean.cred = TRUE, mean.cred.col = "green", xlab)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{over.stat}{
%%     ~~Describe \code{over.stat} here~~
}
  \item{nbreaks}{
%%     ~~Describe \code{nbreaks} here~~
}
  \item{equal.axis}{
%%     ~~Describe \code{equal.axis} here~~
}
  \item{species.names}{
%%     ~~Describe \code{species.names} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{mean.cred}{
%%     ~~Describe \code{mean.cred} here~~
}
  \item{mean.cred.col}{
%%     ~~Describe \code{mean.cred.col} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
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
function (over.stat, nbreaks = 50, equal.axis = FALSE, species.names, 
    col, mean.cred = TRUE, mean.cred.col = "green", xlab) 
{
    if (length(dim(over.stat)) != 3 || dim(over.stat)[1] != dim(over.stat)[2]) 
        stop("Incorrect specification of over.stat.")
    nspec <- dim(over.stat)[1]
    if (missing(species.names)) 
        species.names <- dimnames(over.stat)[[1]]
    over.hist <- apply(over.stat, 1:2, function(x) {
        if (any(is.na(x))) 
            return(NULL)
        else {
            tmp <- hist(x * 100, breaks = nbreaks, plot = FALSE)
            tmp$density <- tmp$density/max(tmp$density)
            tmp$counts <- tmp$counts/max(tmp$counts)
        }
        tmp
    })
    xlim <- lapply(over.hist, function(h) {
        if (is.null(h)) 
            tmp <- matrix(NA, 2, 2)
        else {
            tmp <- cbind(range(h$breaks), range(h$density))
        }
        tmp
    })
    xlim <- matrix(xlim, nspec, nspec)
    if (equal.axis) {
        for (ii in 1:nspec) {
            xlim[, ii] <- rep(list(cbind(range(sapply(xlim[, 
                ii], function(ll) ll[, 1]), na.rm = TRUE), range(sapply(xlim[, 
                ii], function(ll) ll[, 2]), na.rm = TRUE))), 
                nspec)
        }
    }
    if (mean.cred) {
        over.mean <- apply(over.stat, 1:2, mean) * 100
        over.cred <- apply(over.stat, 1:2, quantile, prob = c(0.025, 
            0.975), na.rm = TRUE) * 100
        over.cred <- array(over.cred, dim = c(2, nspec, nspec))
    }
    par(mfcol = c(nspec, nspec), mar = c(1.5, rep(0.5, 3)), oma = rep(4, 
        4))
    for (ci in 1:nspec) {
        for (ri in 1:nspec) {
            plot.new()
            if (ri != ci) {
                plot.window(xlim[ri, ci][[1]][, 1], xlim[ri, 
                  ci][[1]][, 2])
                if (equal.axis) {
                  tmp <- hist(over.stat[ri, ci, ] * 100, breaks = seq(xlim[ri, 
                    ci][[1]][1, 1], xlim[ri, ci][[1]][2, 1], 
                    len = nbreaks + 1), plot = FALSE)
                  tmp$density <- tmp$density/max(tmp$density)
                  over.hist[[ri, ci]] <- tmp
                }
                plot(over.hist[[ri, ci]], add = TRUE, freq = FALSE, 
                  col = col[ci], border = "white")
                if (mean.cred) {
                  abline(v = c(over.mean[ri, ci], over.cred[, 
                    ri, ci]), col = mean.cred.col, lty = c(1, 
                    2, 2), lwd = 2)
                }
            }
            else {
                plot.window(xlim = c(0, 1), ylim = c(0, 1))
            }
            if (ri == 1 && ci == 1) {
                text(x = 0.5, y = 0.5, labels = expression("Niche Overlap: " * 
                  bgroup("(", atop("Row", "Column"), ")")), adj = c(0.5, 
                  0.5), cex = 1)
            }
            box()
            if (ci != ri) 
                axis(side = 1)
            if (ci > 1) 
                axis(side = 2, labels = FALSE)
            if (ci < nspec) 
                axis(side = 4, labels = FALSE)
            if (ri == 1) 
                mtext(text = species.names[ci], side = 3, line = 1, 
                  col = col[ci])
            if (ci == 1) 
                mtext(text = species.names[ri], side = 2, line = 1)
            if (mean.cred && ri == nspec && ci == nspec) {
                legend(x = "center", legend = c("Mean", "95\% Credible Interval "), 
                  lty = c(1, 2), lwd = 2, col = mean.cred.col)
            }
        }
    }
    if (!missing(xlab)) {
        mtext(text = xlab, side = 1, line = 1.5, cex = 1, outer = TRUE)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
