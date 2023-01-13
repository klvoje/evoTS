#' @title Fit multivariate Unbiased Random Walk models to multivariate evolutionary sequence (time-series) data.
#'
#' @description Function to find maximum likelihood solutions to a multivariate Unbiased Random Walk model.
#'
#' @param yy a multivariate evoTS object.
#'
#' @param R the drift matrix. The options are "diagonal" and "symmetric"
#'
#' @param r parameter describing the exponential increase/decrease in rate across time. The options are "fixed", "accel" and "decel".
#'
#' @param method optimization method, passed to function optim. Default is "L-BFGS-B".
#'
#' @param hess logical, indicating whether to calculate standard errors from the Hessian matrix.
#' 
#' @param pool indicating whether to pool variances across samples
#'
#' @param trace logical, indicating whether information on the progress of the optimization is printed.
#'
#' @param iterations the number of times the optimization method is run from different starting points. Default is NULL, meaning the optimization is run once.
#'
#' @param iter.sd defines the standard deviation of the Gaussian distribution from which starting values for the optimization routine is run. Default is 1.
#'
#' @details The function allows the users to test six variants of multivariate Unbiased Random Walk models. There are two options for the structure of the R matrix. A "diagonal" R matrix means the stochastic changes in the traits are assumed to be uncorrelated. A "symmetric" R matrix means the stochastic changes in the traits are assumed to be correlated, i.e. that they are non-independent.
#'
#' There are three options for the 'r' parameter. The "fixed" option means there is no change in the rate of change across time (r = 0). Setting r to "fixed" therefore fits a regular multivariate Unbiased Random Walk. The "decel" and "accel" options make the rate of change (the R matrix) decay (r < 0) and increase (r > 0) exponentially through time, respectively.
#'
#' The function searches - using an optimization routine - for the maximum-likelihood solution for the chosen multivariate Unbiased Random Walk model. The argument 'method' is passed to the 'optim' function and is included for the convenience of users to better control the optimization routine. The the default method (L-BFGS-B) seems to work for most evolutionary sequences.
#'
#' Initial estimates to start the optimization come from maximum-likelihood estimates of the univariate Unbiased Random Walk model (from the paleoTS package) fitted to each time-series separately.
#'
#' It is good practice to repeat any numerical optimization procedure from different starting points. This is especially important for complex models as the log-likelihood surface might contain more than one peak. The number of iterations is controlled by the argument 'iterations'. The function will report the model parameters from the iteration with the highest log-likelihood.
#'
#'@return First part of the output reports the log-likelihood of the model and its AICc score. The second part of the output is the maximum log-likelihood model parameters (ancestral.values, R). The last part of the output gives information about the number of parameters in the model (K), number of samples in the data (n) and number of times the optimization routine was run (iter).
#'
#'@note The models have been implemented to be compatible with the joint parameterization routine in the package paleoTS. The optimization is therefore fit using the actual sample values, with the autocorrelation among samples accounted for in the log-likelihood function. The joint distribution of sample means is multivariate normal, with means and variance-covariances determined by evolutionary parameters and sampling errors.
#'
#'@author Kjetil Lysne Voje
#'
#'@references Revell, L. J. & Harmon, L. Testing quantitative genetic hypotheses about the evolutionary rate matrix for continuous characters. \emph{Evolutionary Ecology Research} 10, 311–331 (2008).
#'@references Voje, K. L. Testing eco‐evolutionary predictions using fossil data: Phyletic evolution following ecological opportunity. \emph{Evolution} 74, 188–200 (2020).
#'
#'@export
#'
#'@examples
#'## Generate an evoTS object by simulating a multivariate dataset
#'x <- sim.multi.URW(30)
#'
#'## Fit a multivariate Unbiased Random Walk model to the data, allowing for correlated changes.
#'fit.multivariate.URW(x, R = "symmetric", r = "fixed")
#'

fit.multivariate.URW<-function(yy, R = "symmetric", r = "fixed", method="L-BFGS-B", hess = FALSE, pool = TRUE, trace=FALSE, iterations=NULL, iter.sd=NULL){

  m <- ncol(yy$xx) # number of traits

  if (R == "symmetric" && r == "fixed") w<-opt.single.R(yy, method = method, hess = hess, pool = pool, trace = trace, iterations = iterations, iter.sd = iter.sd)

  if (R == "diag" && r == "fixed") w<-opt.single.R.zero.corr(yy, method = method, hess = hess, pool = pool, trace = trace, iterations = iterations, iter.sd = iter.sd)

  if (R == "symmetric" && r == "accel") w<-opt.accel.single.R(yy, method = method, hess = hess, pool = pool, trace = trace,  iterations = iterations, iter.sd = iter.sd)

  if (R == "symmetric" && r == "decel") w<-opt.decel.single.R(yy, method = method, hess = hess, pool = pool, trace = trace, iterations = iterations, iter.sd = iter.sd)

  if (R == "diag" && r == "accel") w<-opt.accel.single.R.zero.corr(yy, method = method, hess = hess, pool = pool, trace = trace, iterations = iterations, iter.sd = iter.sd)

  if (R == "diag" && r == "decel") w<-opt.decel.single.R.zero.corr(yy, method = method, hess = hess, pool = pool, trace = trace, iterations = iterations, iter.sd = iter.sd)

  return(w)
}
