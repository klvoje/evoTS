## Update of the evoTS package
 
## summary of the new changes

#Bug fixes
- fixed a bug in functions running multivariate models with iterations: The functions are now able to provide output from the best model even if the initial parameters failed to work for one or more of the iterations. 

# Other changes
- Added functionality for running iterations of the multivariate unbiased random walk moded with a shift point (i.e., where different sections of a time-series are described by separate R matrices).

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

