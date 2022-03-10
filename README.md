
<!-- README.md is generated from README.Rmd. Please edit that file -->

# evoTS - analysis of evolutionary sequences of phenotypic change

<!-- badges: start -->
<!-- badges: end -->

The `evoTS` package facilitates univariate and multivariate analysis of
evolutionary sequences of phenotypic change. The goal of the `evoTS`
package is to offer a large range of evolutionary models to enable more
detailed studies of evolutionary changes within lineages.

The `evoTS` package extends the modeling framework available in the
<a href="https://CRAN.R-project.org/package=paleoTS"> `paleoTS`
package</a>. All model-fitting procedures in `evoTS` have been
implemented to mirror the user experience from `paleoTS`. For example,
all univariate models implemented in `evoTS` can be fitted to a
`paleoTS` object, i.e. the data format used in `paleoTS`. The fit of all
univariate models available in `paleoTS` and `evoTS` are directly
comparable using the reported AICc.

`evoTS` contains functions that allow for fitting different models to
separate parts of an evolutionary sequence. Functions for investigating
likelihood surfaces of fitted models are also included.

`evoTS` contains a range of multivariate models, including different
versions of multivariate unbiased random walks and Ornstein-Uhlenbeck
processes. Together, these models allow the user to test various
hypotheses of adaptation using phenotypic time-series.

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
