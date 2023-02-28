# evoTS 1.0.1 -> Current version on CRAN

## Bug fixes

- Fixed a bug in functions running multivariate models with iterations: The functions are now providing output from the best model even if the initial parameters failed to work for one or more of the iterations. 

- Fixed a bug in how the variance-covariance matrix of the decelerated and accelerated models were defined.    

## Other changes

- Added functionality for running iterations of the multivariate unbiased random walk model with a shift point (i.e., where different sections of a time-series are described by separate R matrices).

## Adjustments currently only implemented the the development version on GitHuB

- Implemented a change in the box constraints in the L-BFGS-B method for the multivariate accelerated and decelerated models. 
- Updated the list of messages to print as part of the output based on how the 'optim' function may fail while searching for the maximum likelihood for multivariate unbiased random walk models.  

