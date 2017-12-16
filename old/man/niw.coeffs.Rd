\name{niw.coeffs}
\alias{niw.coeffs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
niw.coeffs(X, lambda, kappa, Psi, nu)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{kappa}{
%%     ~~Describe \code{kappa} here~~
}
  \item{Psi}{
%%     ~~Describe \code{Psi} here~~
}
  \item{nu}{
%%     ~~Describe \code{nu} here~~
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
function (X, lambda, kappa, Psi, nu) 
{
    d <- ncol(X)
    N <- nrow(X)
    if (missing(kappa)) 
        kappa <- 0
    if (missing(Psi)) 
        Psi <- 0
    if (missing(nu)) 
        nu <- ncol(X) + 1
    Xbar <- colMeans(X)
    S <- t(X) - Xbar
    S <- S \%*\% t(S)
    Psi2 <- Psi + S
    lambda2 <- N * Xbar
    if (kappa != 0) {
        Psi2 <- Psi2 + (N * kappa)/(N + kappa) * (Xbar - lambda) \%*\% 
            t(Xbar - lambda)
        lambda2 <- lambda2 + kappa * lambda
    }
    lambda2 <- lambda2/(N + kappa)
    nu2 <- N + nu - (kappa == 0)
    kappa2 <- N + kappa
    list(lambda = lambda2, kappa = kappa2, Psi = Psi2, nu = nu2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
