
<!-- README.md is generated from README.Rmd. Please edit that file -->

# evoTS

<!-- badges: start -->
<!-- badges: end -->

The evoTS package facilitates analysis of evolutionary sequences of
phenotypic change. The package extends the modeling framework in the
paleoTS package (Hunt 2006; 2008a; 2008b; Hunt et al. 2008; 2010; 2015)
and contains functions to fit a large range of univariate and
multivariate evolutionary models not implemented in paleoTS. The goal of
the evoTS package is to offer a larger range of evolutionary models to
te study of evolutionary changes within lineages.

All model-fitting procedures in evoTS have been implemented to mirror
the user experience from paleoTS in order to make the package easy to
use for users with experince from paleoTS. For example, all the
univariate models implemented in evoTS can be fitted to a paleoTS object
(i.e. the data format used in paleoTS) and the fit of all univariate
models in paleoTSand evoTS are directly comparable using the reported
AICc. evoTS also contain functions that allow for more flexibility in
fitting combinations of univariate models to different parts of a
time-series. Functions for investigating likelihood surfaces of fitted
models have also been developed.

evoTS also contains a range of multivariate models that can be fitted to
evolutionary sequences (time-series). These models include different
versions of multivariate Unbiased Random Walks and Ornsten-Uhlenbeck
processes.

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
