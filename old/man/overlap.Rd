\name{overlap}
\alias{overlap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
overlap(niche.par, nreps, nprob, alpha = 0.95, species.names, norm.redraw = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{niche.par}{
%%     ~~Describe \code{niche.par} here~~
}
  \item{nreps}{
%%     ~~Describe \code{nreps} here~~
}
  \item{nprob}{
%%     ~~Describe \code{nprob} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{species.names}{
%%     ~~Describe \code{species.names} here~~
}
  \item{norm.redraw}{
%%     ~~Describe \code{norm.redraw} here~~
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
function (niche.par, nreps, nprob, alpha = 0.95, species.names, 
    norm.redraw = TRUE) 
{
    niso <- ncol(niche.par[[1]]$mu)
    nspec <- length(niche.par)
    nlevels <- length(alpha)
    if (missing(species.names)) 
        species.names <- names(niche.par)
    mu <- matrix(NA, niso, nspec)
    x <- array(NA, dim = c(niso, nprob, nspec))
    C <- array(NA, dim = c(niso, niso, nspec))
    Zsq <- array(NA, dim = c(nprob, nspec - 1, nlevels))
    qlevels <- qchisq(alpha, df = niso)
    qlevels <- array(rep(qlevels, each = nprob * (nspec - 1)), 
        dim = c(nprob, nspec - 1, nlevels))
    if (!norm.redraw) 
        z <- matrix(rnorm(niso * nprob), niso, nprob)
    over <- array(NA, dim = c(nspec, nspec, nreps, nlevels))
    dimnames(over) <- list(species.names, species.names, NULL, 
        paste0(alpha * 100, "\%"))
    names(dimnames(over)) <- c("Species A", "Species B", "", 
        "alpha")
    ind <- sapply(niche.par, function(mv) {
        sample(nrow(mv$mu), nreps, replace = TRUE)
    })
    for (jj in 1:nreps) {
        for (ii in 1:nspec) {
            mu[, ii] <- niche.par[[ii]]$mu[ind[jj, ii], ]
            C[, , ii] <- t(chol(niche.par[[ii]]$Sigma[, , ind[jj, 
                ii]]))
            if (norm.redraw) 
                z <- matrix(rnorm(niso * nprob), niso, nprob)
            x[, , ii] <- C[, , ii] \%*\% z + mu[, ii]
        }
        for (ii in 1:nspec) {
            Z <- backsolve(C[, , ii], matrix(x[, , -ii], nrow = niso) - 
                mu[, ii], upper.tri = FALSE)
            Zsq[] <- colSums(Z^2)
            over[-ii, ii, jj, ] <- colMeans(Zsq < qlevels)
        }
    }
    if (dim(over)[4] == 1) 
        over <- over[, , , 1]
    over
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
