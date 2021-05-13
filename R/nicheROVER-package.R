#' (Niche) (R)egion and Niche (Over)lap Metrics for Multidimensional Ecological Niches.
#'
#' This package uses a probabilistic method to calculate niche regions and pairwise niche overlap using multidimensional niche indicator data (e.g., stable isotopes, environmental variables, etc.). The niche region is defined as the joint probability density function of the multidimensional niche indicators at a user-defined probability alpha (e.g., 95%).  Uncertainty is accounted for in a Bayesian framework, and the method can be extended to three or more indicator dimensions.  It provides directional estimates of niche overlap, accounts for species-specific distributions in multivariate niche space, and produces unique and consistent bivariate projections of the multivariate niche region. See Swanson et al. (2014) for a detailed description and worked example below using fish stable isotope data.
#'
#' @references Swanson, H.K., Lysy, M., Stasko, A.D., Power, M., Johnson, J.D., and Reist, J.D. "A new probabilistic method for quantifying n-dimensional ecological niches and niche overlap." *Ecology: Statistical Reports* 96:2 (2015): 318-324. <https://www.ncbi.nlm.nih.gov/pubmed/26240852>.
#'
#' @example examples/nicheROVER-package.R
#'
#' @docType package
#' @importFrom graphics abline axis box hist legend lines mtext par plot plot.new plot.window points segments text
#' @importFrom stats density pbeta qchisq quantile rbeta rchisq rnorm
#' @importFrom utils combn
#' @name nicheROVER
NULL
