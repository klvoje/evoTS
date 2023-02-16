# evoTS 1.0.1

## Bug fixes

- Fixed a bug in functions running multivariate models with iterations: The functions are now  providing output from the best model even if the initial parameters failed to work for one or more of the iterations. 

## Other changes

- Added functionality for running iterations of the multivariate unbiased random walk model with a shift point (i.e., where different sections of a time-series are described by separate R matrices).