# evoTS 1.0.3

## Bug fixes

- Fixed a bug in how the variance-covariance matrix of the decelerated and accelerated models were defined. 

## Other changes

- Implemented a change in the box constraints in the L-BFGS-B method for the multivariate accelerated and decelerated models. 
- Updated the list of messages to print as part of the output based on how the 'optim' function may fail while searching for the maximum likelihood for multivariate unbiased random walk models.
- The output format for univariate models in evoTS has been updated to ensure compatibility with the latest version (0.6.1) of paleoTS, which has altered how it displays such output.

## Adjustments currently only implemented the the development version on GitHuB

The version of evoTS on GitHub is currently identical to the latest (v.1.0.3), the version on CRAN.

