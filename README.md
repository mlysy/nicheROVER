# nicheROVER: Niche Region and Niche Overlap Metrics for Multidimensional Ecological Niches

*Martin Lysy, Ashley Stasko, Heidi Swanson*

---

### Description

This **R** package uses a probabilistic method to calculate niche regions and pairwise niche overlap using multidimensional niche indicator data (e.g., stable isotopes, environmental variables, etc.). The niche region is defined as the joint probability density function of the multidimensional niche indicators at a user-defined probability alpha (e.g., 95%).  Uncertainty is accounted for in a Bayesian framework, and the method can be extended to three or more indicator dimensions.  It provides directional estimates of niche overlap, accounts for species-specific distributions in multivariate niche space, and produces unique and consistent bivariate projections of the multivariate niche region.  The article by [Swanson et al. (Ecology, 2015)](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/14-0235.1) provides a detailed description of the methodology.  See the [package vignette](http://htmlpreview.github.com/?https://github.com/mlysy/nicheROVER/master/inst/doc/ecol-vignette.html) for a worked example using fish stable isotope data.

### Installation

The CRAN version is available [here](https://CRAN.R-project.org/package=nicheROVER).

For the more recent GitHub version (in particular, including niche size calculations), first install the **R** package [**`devtools`**](https://CRAN.R-project.org/package=devtools), then run
```{r}
devtools::install_github("mlysy/nicheROVER", ref = "master")
```

### Usage

A step-by-step walkthrough of a basic analysis is provided in the package [vignette](http://htmlpreview.github.com/?https://github.com/mlysy/nicheROVER/master/doc/ecol-vignette.html).
