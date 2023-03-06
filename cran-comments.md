## Update of the evoTS package
 
## summary of the new changes

#Bug fixes
- Fixed a bug in how the variance-covariance matrix of the decelerated and accelerated models was defined. 

# Other changes
- Implemented a change in the box constraints in the L-BFGS-B method for the multivariate accelerated and decelerated models. 
- Updated the list of messages to print as part of the output based on how the 'optim' function may fail while searching for the maximum likelihood for multivariate unbiased random walk models. 

## Test environments
Windows Server 2022, R-devel, 64 bit
Fedora Linux, R-devel, clang, gfortran
Ubuntu Linux 20.04.1 LTS, R-release, GCC
Development version of R
Release version of R

## R CMD check results
0 errors v | 0 warnings v | 0 notes v

## Downstream dependencies
There are currently no downstream dependencies for this package

