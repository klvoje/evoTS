#' @title Fit separate multivariate Unbiased Random Walk models to two different segments of a multivariate evolutionary sequence (time-series).
#'
#' @description Function to find maximum likelihood solutions for multivariate Unbiased Random Walk models fitted to two different segments of a multivariate evolutionary sequence (time-series).
#'
#' @param yy a multivariate evoTS object.
#'
#' @param minb the minimum number of samples within a segment to consider.
#'
#' @param hess logical, indicating whether to calculate standard errors from the Hessian matrix.
#' 
#' @param pool indicating whether to pool variances across samples
#'
#' @param shift.point the sample in the time series that represents the first sample in the second segment.
#'
#' @param method optimization method, passed to function optim. Default is "L-BFGS-B".
#'
#' @param trace logical, indicating whether information on the progress of the optimization is printed.
#'
#' @param iterations the number of times the optimization method is run from different starting points. Default is NULL, meaning the optimization is run once.
#'
#' @param iter.sd defines the standard deviation of the Gaussian distribution from which starting values for the optimization routine is run. Default is 1.
#'
#' @details The function searches - using an optimization routine - for the maximum-likelihood solution for a multivariate Unbiased Random Walk model ti two non-overlapping segments in the time series.
#'
#' The argument 'method' is passed to the 'optim' function and is included for the convenience of users to better control the optimization routine. The the default method (L-BFGS-B) seems to work for most evolutionary sequences.
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
#'
#'@export
#'
#'@examples
#'## Generate an evoTS object by simulating a multivariate dataset
#'x <- sim.multi.URW(60)
#'
#'## Fit two multivariate Unbiased Random Walk models to separate parts of the time-series.
#'fit.multivariate.URW.shift(x, shift.point = 31)

fit.multivariate.URW.shift<-function (yy, minb = 10, hess = FALSE, pool = TRUE, shift.point = NULL, method = "L-BFGS-B", trace=FALSE, iterations=NULL, iter.sd=NULL)
{

  ng<-2 # The program is currently constrained to fitting two R matrices.
  ns <- nrow(yy$xx)
  m <- ncol(yy$xx) # number of traits
  
  if (pool==TRUE) { 
    for (i in 1:m){
      
      tmp<-paleoTS::as.paleoTS(yy$xx[,i], yy$vv[,i], yy$nn[,i], yy$tt[,i])
      tmp<- paleoTS::pool.var(tmp, ret.paleoTS = TRUE)
      yy$vv[,i]<-tmp$vv
    }
  }

  if(is.numeric(shift.point) == TRUE) GG <-shift.point else GG <- shifts(ns, ng, minb = minb)
  GG<-as.matrix(GG)

  #Define number of shift points:
  nc <- ncol(GG)
    cat("Total # hypotheses: ", nc, "\n")
  #Create empty list
  wl <- list()
  #create array with length = to shift points and where every entry = -Inf
  logl <- array(-Inf, dim = nc)

  #start loop for estimating maximum likelihood parameters for each data set defined by the shift points
  for (i in 1:nc) {
      cat(i, " ")
    #defines which data point in the time series that belong to each of the two sets
    gg <- shift2gg(GG[, i], ns)
        w <- opt.multi.R(yy, gg, hess = hess, pool = pool, method = method, trace = trace, iterations=iterations, iter.sd=iter.sd)
    logl[i] <- w$logL
    wl[[i]] <- w
  }

    cat("\n")
  winner <- which.max(logl)
  ww <- wl[[winner]]
  ss <- GG[, winner]
  names(ss) <- paste("shift", 1:(ng - 1), sep = "")
  ww$parameters <- append(ww$parameters, ss)
  ww$all.logl <- logl
  ww$GG <- GG
  return(ww)
}
