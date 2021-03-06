% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/niche.size.R
\name{niche.size}
\alias{niche.size}
\title{Calculate the size of an elliptical niche region.}
\usage{
niche.size(Sigma, alpha = 0.95)
}
\arguments{
\item{Sigma}{Variance matrix for normally distributed niche axes.}

\item{alpha}{Probabilistic niche size.}
}
\value{
Hypervolume niche size.
}
\description{
Calculate the size of an elliptical niche region.
}
\details{
For a given niche region \eqn{N_R}, the niche size is defined as the hypervolume of this region: \eqn{N_S = \int_{x \in N_R} d x}.
}
\examples{
# for each species, size of 95\% niche region using sample variance
tapply(1:nrow(fish), fish$species, function(ind) {
  X <- fish[ind,2:4] # all measurements for given species
  Sighat <- var(X) # sample variance
  niche.size(Sigma = Sighat, alpha = .95)
})
}
