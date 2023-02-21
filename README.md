
<!-- README.md is generated from README.Rmd. Please edit that file -->

# evoTS - Analyses of evolutionary time-series <img src="vignettes/evoTS.png" width="120" align="right" />

<!-- badges: start -->

[![](https://www.r-pkg.org/badges/version/evoTS?color=blue)](https://cran.r-project.org/package=evoTS)
[![](http://cranlogs.r-pkg.org/badges/grand-total/evoTS?color=blue)](https://cran.r-project.org/package=evoTS)
<!-- badges: end -->

<https://cran.r-project.org/package=evoTS>

The `evoTS` package facilitates univariate and multivariate analyses of
phenotypic change within lineages.

The `evoTS` package extends the modeling framework available in the
<a href="https://CRAN.R-project.org/package=paleoTS"> `paleoTS`
package</a>. All model-fitting procedures in `evoTS` have been
implemented to mirror the user experience from `paleoTS`. For example,
all univariate models implemented in `evoTS` can be fitted to a
`paleoTS` object, i.e. the data format used in `paleoTS`. The fit of all
univariate models available in `paleoTS` and `evoTS` are directly
comparable using the reported AICc.

`evoTS` contains functions that allow for fitting different models to
separate parts of an evolutionary sequence (mode-shift models).
Functions for investigating likelihood surfaces of fitted models are
also included.

Multivariate models implemented in `evoTS` include different versions of
multivariate unbiased random walks and Ornstein-Uhlenbeck processes.
These multivariate models allow the user to test a variety of hypotheses
of adaptation and evolution using phenotypic time-series.

## The development version

The GitHub repository contains a copy of the current development version
of the R package evoTS. This version is as recent as or more recent than
the official release of evoTS on the Comprehensive R Archive Network
(CRAN), which is available at
<https://cran.r-project.org/package=evoTS>.

## Where is the official (stable) release?

For the most recent official and stable release of `evoTS`, see
<https://cran.r-project.org/package=evoTS>

## Installation

### Installing the official release

``` r
## Installing from CRAN
> install.packages("evoTS")

> library(evoTS)
```

### Installing the the development version from GitHub

``` r
## Installing from GitHub
> install.packages("devtools")

> devtools::install_github("klvoje/evoTS")

> library(evoTS)
```

## Documentation

The <a href="https://klvoje.github.io/evoTS/index.html">package
website</a> contains a vignette (detailed walk-through) on how to use
the various features of the `evoTS` package.
