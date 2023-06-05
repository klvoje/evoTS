# evoTS 1.0.2

## Bug fixes

- Fixed a bug in how the variance-covariance matrix of the decelerated and accelerated models were defined. 

## Other changes

- Implemented a change in the box constraints in the L-BFGS-B method for the multivariate accelerated and decelerated models. 
- Updated the list of messages to print as part of the output based on how the 'optim' function may fail while searching for the maximum likelihood for multivariate unbiased random walk models. 

## Adjustments currently only implemented the the development version on GitHuB

- Updated the list of messages to print as part of the output based on how the 'optim' function may fail while searching for the maximum likelihood for multivariate unbiased random walk models with a shift.
- Fixed a bug in how the maximum likelihood parameters were reported in the multivariate unbiased random walk model with a shift. The bug only affected the output when a shift point was not defined while iterations was defined. 
- Changed some of the initial starting values in one of the functions running multivariate OU models (fit.multivariate.OU) to make all starting values identical across all functions and iteration options. The starting values of off-diagonal elements in A and R are now always 0. 
- Changed the function fit.all.univariate to prevent it from failing when one of the variances is 0 and the sample size is 1.

