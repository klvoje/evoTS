
<!-- README.md is generated from README.Rmd. Please edit that file -->

# evoTS - univariate and multivariate models of phenotypic evolution within lineages

<!-- badges: start -->
<!-- badges: end -->

The `evoTS` package facilitates analysis of evolutionary sequences of
phenotypic change. The package extends the modeling framework in the
<a href="https://cran.r-project.org/web/packages/paleoTS/index.html">
`paleoTS` package</a> and contains functions to fit univariate and
multivariate evolutionary models not implemented in `paleoTS.` The goal
of the `evoTS` package is to offer a larger range of evolutionary models
to the study of evolutionary changes within lineages.

All model-fitting procedures in `evoTS` have been implemented to mirror
the user experience from `paleoTS`. All the univariate models
implemented in `evoTS` can be fitted to a `paleoTS` object (i.e. the
data format used in `paleoTS`). The fit of all univariate models
available in `paleoT` Sand `evoTS` are directly comparable using the
reported AICc.

`evoTS` contains functions that allow for fitting different models to
separate parts of an evolutionary sequence. Functions for investigating
likelihood surfaces of fitted models are also included.

`evoTS` contains a range of multivariate models that can be fitted to
evolutionary sequences (time-series), including different versions of
multivariate unbiased random walks and Ornsten-Uhlenbeck processes.
Together, these models allow the user to test various hypotheses of
adaptation using phenotypic time-series

## Installation

The evoTS package is available on GitHub and can be installed using
devtools:

``` r
install.packages("devtools")
devtools::install_github("klvoje/evoTS")
```

## Documentation

The <a href="https://klvoje.github.io/evoTS/">package website</a>
contains a vignette (detailed walk-through) on how to use the various
features of the `evoTS` package.
