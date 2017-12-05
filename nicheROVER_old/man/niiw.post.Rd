\name{niiw.post}
\alias{niiw.post}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
niiw.post(nsamples, X, lambda, Omega, Psi, nu, mu0 = lambda, burn)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nsamples}{
%%     ~~Describe \code{nsamples} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{Omega}{
%%     ~~Describe \code{Omega} here~~
}
  \item{Psi}{
%%     ~~Describe \code{Psi} here~~
}
  \item{nu}{
%%     ~~Describe \code{nu} here~~
}
  \item{mu0}{
%%     ~~Describe \code{mu0} here~~
}
  \item{burn}{
%%     ~~Describe \code{burn} here~~
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
function (nsamples, X, lambda, Omega, Psi, nu, mu0 = lambda, 
    burn) 
{
    d <- ncol(X)
    N <- nrow(X)
    Xbar <- colMeans(X)
    S <- t(X) - Xbar
    S <- S \%*\% t(S)
    mu <- rep(0, d)
    Sigma <- matrix(0, d, d)
    if (missing(burn)) 
        burn <- 0.1
    if (burn < 1) 
        burn <- floor(nsamples * burn)
    mu <- mu0
    mu.out <- matrix(NA, nsamples, d)
    Sigma.out <- array(NA, dim = c(d, d, nsamples))
    for (ii in (-burn + 1):nsamples) {
        Psi2 <- N * ((mu - Xbar) \%o\% (mu - Xbar)) + S + Psi
        nu2 <- N + nu
        Sigma <- matrix(rwish(1, Psi2, nu2, inv = TRUE), d, d)
        Sigma2 <- Sigma/N
        if (!all(Omega == 0)) {
            B <- Sigma2 \%*\% solve(Omega + Sigma2)
        }
        else B <- matrix(0, d, d)
        IB <- diag(d) - B
        lambda2 <- c(IB \%*\% Xbar)
        if (!all(Omega == 0)) 
            lambda2 <- lambda2 + c(B \%*\% lambda)
        Omega2 <- IB \%*\% Sigma2
        mu <- c(rmvnorm(1, lambda2, Omega2))
        if (ii > 0) {
            mu.out[ii, ] <- mu
            Sigma.out[, , ii] <- Sigma
        }
    }
    list(mu = mu.out, Sigma = Sigma.out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
