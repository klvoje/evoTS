## Update of the evoTS package
 
## summary of the new changes

#Bug fixes
- Fixed a bug in how the variance-covariance matrix of the decelerated and accelerated models was defined. 

# Other changes
- Implemented a change in the box constraints in the L-BFGS-B method for the multivariate accelerated and decelerated models. 
- Updated the list of messages to print as part of the output based on how the 'optim' function may fail while searching for the maximum likelihood for multivariate unbiased random walk models. 
- The output format for univariate models in evoTS has been updated to ensure compatibility with the latest version (0.6.1) of paleoTS, which has altered how it displays such output.  


## Test environments
Release version of R
Development version of R
windows R-* (any version) windows-latest on GitHub
atlas, R-devel (2024-06-18 r86781), Fedora Linux 38 (Container Image)

## R CMD check results
0 errors v | 0 warnings v | 0 notes v

## Downstream dependencies
There are currently no downstream dependencies for this package

