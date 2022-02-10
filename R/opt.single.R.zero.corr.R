#' @title Fit multivariate Unbiased Random Walk model with uncorrelated trait changes.
#'
#' @description Function to find maximum likelihood solution to a multivariate Unbiased Random Walk model with uncorrelated trait changes.
#'
#' @param yy a multivariate evoTS object.
#'
#' @param method optimization method, passed to function optim. Default is "L-BFGS-B".
#'
#' @param hess logical, indicating whether to calculate standard errors from the Hessian matrix.
#'
#' @param trace logical, indicating whether information on the progress of the optimization is printed.
#'
#' @param iterations the number of times the optimization method is run from different starting points. Default is NULL, meaning the optimization is run once.
#'
#' @param iter.sd defines the standard deviation of the Gaussian distribution from which starting values for the optimization routine is run. Default is 1.
#'
#' @details The function searches - using an optimization routine - for the maximum-likelihood solution for a multivariate Unbiased Random Walk model with uncorrelated trait changes.
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
#'## Generate an evoTS objects by simulating a multivariate dataset
#'x <- sim.multi.URW(30)
#'
#'## Fit a multivariate Unbiased Random Walk model with uncorrelated trait changes.
#'opt.single.R.zero.corr(x)

opt.single.R.zero.corr<-function (yy, method = "L-BFGS-B", hess = FALSE, trace=FALSE, iterations=NULL, iter.sd=NULL)
{

  n <- nrow(yy$xx) #number of samples/populations
  m <- ncol(yy$xx) # number of traits

  C<- outer(yy$tt[,1], yy$tt[,1], FUN = pmin)  #Create distance matrix

  X <- yy$xx # Character matrix with dimensions n * m
  y <- as.matrix(as.vector(X)) # Vectorized version of X

  # Define initial parameter values for the optimization routine
  init.trait.var<-apply(yy$xx,2,var)
  anc.values<-yy$xx[1,]

  init.par<-c(init.trait.var, anc.values)
  lower.limit<-c(rep(0,length(init.trait.var)), rep(NA, length(anc.values)))

  ### Start iterations from different starting values
  if (is.numeric(iterations)) {
    if(is.numeric(iter.sd) == FALSE) iter.sd <-1
    #if(is.numeric(max.attemps) == FALSE) max.attemps <-100000
    log.lik.tmp<-rep(NA, iterations)
    www<-list()

    for (k in 1:iterations){

      init.par_temp<-init.par
      init.par<-rnorm(length(init.par_temp), init.par_temp, iter.sd)


      if (method == "L-BFGS-B")  {
        www[[k]]<-optim(init.par, fn = logL.joint.single.R.zero.corr, C = C, y = y, m = m, n = n, anc.values = anc.values, yy = yy,
                        control = list(fnscale = -1, maxit=10000, trace = trace), method = "L-BFGS-B", hessian = hess, lower = lower.limit)
      }

      if (method == "Nelder-Mead")  {
        www[[k]]<-optim(init.par, fn = logL.joint.single.R.zero.corr, C = C, y = y, m = m, n = n, anc.values = anc.values, yy = yy,
                        control = list(fnscale = -1, maxit=10000, trace = trace), method = "Nelder-Mead" , hessian = hess)
      }
      if (method == "SANN")  {
        www[[k]]<-optim(init.par, fn = logL.joint.single.R.zero.corr, C = C, y = y, m = m, n = n, anc.values = anc.values, yy = yy,
                        control = list(fnscale = -1, maxit=10000, trace = trace), method = "SANN" , hessian = hess, lower = lower.limit)
      }
      log.lik.tmp[k]<-www[[k]]$value

    }
    for (j in 1:iterations){
      if(max(na.exclude(log.lik.tmp)) == www[[j]]$value) best.run<-www[[j]]
    }
  }


  ##### Start of non-iteration routine #####

  if (is.numeric(iterations) == FALSE) {


    if (method == "L-BFGS-B")  {
      w<-optim(init.par, fn = logL.joint.single.R.zero.corr, C = C, y = y, m = m, n = n, anc.values = anc.values, yy = yy,
               control = list(fnscale = -1, maxit=10000, trace = trace), method = "L-BFGS-B", hessian = hess, lower = lower.limit)
    }

    if (method == "Nelder-Mead")  {
      w<-optim(init.par, fn = logL.joint.single.R.zero.corr, C = C, y = y, m = m, n = n, anc.values = anc.values, yy = yy,
               control = list(fnscale = -1, maxit=10000, trace = trace), method = "Nelder-Mead" , hessian = hess)
    }
    if (method == "SANN")  {
      w<-optim(init.par, fn = logL.joint.single.R.zero.corr, C = C, y = y, m = m, n = n, anc.values = anc.values, yy = yy,
               control = list(fnscale = -1, maxit=10000, trace = trace), method = "SANN" , hessian = hess, lower = lower.limit)
    }


    if (w$convergence == 1) print("Model did not converge.")
    if (w$convergence == 0) print("Model converged successfully.")
  }
  # number of parameters
  K <- length(init.par) #parameters in the R matrices + ancestral values for each trait


  if (is.numeric(iterations) ==TRUE) {
     iter<-NA
    w<-best.run
  }

  if (hess) {
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
    SE_R<-matrix(0, nrow=m, ncol=m)
    diag(SE_R)<-w$se[1:m]
    SE.R<-SE_R
  }
  if (hess == FALSE) {
    SE.R<-NA
  }

  R<-matrix(0, nrow=m, ncol=m)
  diag(R)<-w$par[1:m]
  R<-t(R)%*%R

  if (is.numeric(iterations) == FALSE) {
    iter<-NA
  }

  ancestral.values<-w$par[(length(init.trait.var)+1) : length(init.par)]


  wc<-as.evoTS.multi.BW.fit(modelName = "Multivariate Random walk (R matrix with zero off-diagonal elements)", logL = w$value, ancestral.values = ancestral.values, R = R, SE.R = SE.R,
                              method = "Joint", K = K, n = length(yy$xx[,1]), iter=iter)
  return(wc)
}