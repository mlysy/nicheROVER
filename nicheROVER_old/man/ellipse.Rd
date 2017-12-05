\name{ellipse}
\alias{ellipse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ellipse(mu, V, alpha = 0.95, n = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{V}{
%%     ~~Describe \code{V} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
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
function (mu, V, alpha = 0.95, n = 100) 
{
    tmp <- eigen(V)
    hlen <- sqrt(qchisq(alpha, df = 2) * tmp$val)
    theta <- atan2(tmp$vec[2, 1], tmp$vec[1, 1])
    t <- seq(0, 2 * pi, len = n + 1)
    x <- hlen[1] * cos(t)
    y <- hlen[2] * sin(t)
    alpha <- atan2(y, x)
    rad <- sqrt(x^2 + y^2)
    cbind(x = rad * cos(alpha + theta) + mu[1], y = rad * sin(alpha + 
        theta) + mu[2])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
