
<!-- README.md is generated from README.Rmd. Please edit that file -->

# evoTS - Analyses of evolutionary time-series

<!-- badges: start -->
<!-- badges: end -->

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

## Installation

The evoTS package is available on GitHub and can be installed using
devtools:

``` r
install.packages("devtools")
devtools::install_github("klvoje/evoTS")
```

## Documentation

The <a href="https://klvoje.github.io/evoTS/index.html">package
website</a> contains a vignette (detailed walk-through) on how to use
the various features of the `evoTS` package.
