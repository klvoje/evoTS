#' @title Fit separate multivariate Unbiased Random Walk models to two different segments of a multivariate evolutionary sequence (time-series).
#'
#' @description This function is used internally by evoTS and will generally not be used directly by users.
#'
#' @param yy a multivariate evoTS object.
#'
#' @param gg numeric vector indicating membership of each sample in the two segments.
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
#' @details This function is used internally by evoTS and will generally not be used directly by users.
#'
#' The function searches - using an optimization routine - for the maximum-likelihood solution for a multivariate Unbiased Random Walk model.
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

opt.multi.R<-function (yy, gg, method, hess = FALSE, pool = TRUE, trace, iterations, iter.sd)
{

  n <- nrow(yy$xx) # number of samples/populations
  m <- ncol(yy$xx) # number of traits

  states <- length(unique(gg))

  #define vectors based on length of first and second segments in the time series
  seg.1 <- which(gg == 1)
  seg.2 <- which(gg == 2)

  #Extract the time vector for each random walks:
  tt.seg.1 <- yy$tt[,1][seg.1]
  tt.seg.2 <- yy$tt[,1][seg.2]

  #Create distance matrices
  C.1<-matrix(0, n, n)
  C.1[seg.1, seg.1]<- outer(tt.seg.1, tt.seg.1, FUN = pmin)
  C.2<-matrix(0, n, n)
  C.2 [seg.2, seg.2] <- outer(tt.seg.2, tt.seg.2, FUN = pmin)

  C <- list(C.1 = C.1,  C.2= C.2)

  X <- yy$xx # Character matrix with dimensions n * m
  y <- as.matrix(as.vector(X)) # Vectorized version of X


  ### Estimate two full R matrices

  # Define initial parameter values for the optimization routine
  init.trait.var<-apply(yy$xx,2,var)
  temp.matrix<-cov(as.matrix(yy$xx))
  init.cov.traits<-unique(temp.matrix[row(temp.matrix)!=col(temp.matrix)])
  anc.values<-yy$xx[1,]

  init.par<-c(init.trait.var, init.trait.var, init.cov.traits, init.cov.traits, anc.values)
  lower.limit<-c(rep(0,length(init.trait.var)), rep(0,length(init.trait.var)),  rep(NA,length(init.cov.traits)), rep(NA, length(init.cov.traits)), rep(NA, length(anc.values)))

  ### Start iterations from different starting values
  if (is.numeric(iterations)) {
    if(is.numeric(iter.sd) == FALSE) iter.sd <-1
    log.lik.tmp<-rep(NA, iterations)
    www<-list()

    for (k in 1:iterations){

      init.par_temp<-init.par
      init.par<-rnorm(length(init.par_temp), init.par_temp, iter.sd)

        if (method == "L-BFGS-B")  {
          www[[k]]<-optim(init.par, fn = logL.joint.multi.R, C = C, y = y, m = m, n = n, anc.values = anc.values, yy = yy,
                   control = list(fnscale = -1, maxit=10000, trace = trace), method = "L-BFGS-B", hessian = hess, lower = lower.limit)
        }

        if (method == "Nelder-Mead")  {
          www[[k]]<-optim(init.par, fn = logL.joint.multi.R, C = C, y = y, m = m, n = n, anc.values = anc.values, yy = yy,
                   control = list(fnscale = -1, maxit=10000, trace = trace), method = "Nelder-Mead" , hessian = hess)
        }
        if (method == "SANN")  {
          www[[k]]<-optim(init.par, fn = logL.joint.multi.R, C = C, y = y, m = m, n = n, anc.values = anc.values, yy = yy,
                   control = list(fnscale = -1, maxit=10000, trace = trace), method = "SANN" , hessian = hess, lower = lower.limit)
        }
      log.lik.tmp[k]<-www[[k]]$value

    }
    for (j in 1:iterations){
      if(max(na.exclude(log.lik.tmp)) == www[[j]]$value) best.run<-www[[j]]
    }
    w<-best.run
  }


  ##### End of iteration routine #####

  ##### Start of non-iteration routine #####

  if (is.numeric(iterations) == FALSE) {


   if (method == "L-BFGS-B")  {
  w<-optim(init.par, fn = logL.joint.multi.R, C = C, y = y, m = m, n = n, anc.values = anc.values,  yy = yy,
             control = list(fnscale = -1, maxit=10000, trace = trace), method = "L-BFGS-B", hessian = hess, lower = lower.limit)
                              }

  if (method == "Nelder-Mead")  {
  w<-optim(init.par, fn = logL.joint.multi.R, C = C, y = y, m = m, n = n, anc.values = anc.values,  yy = yy,
             control = list(fnscale = -1, maxit=10000, trace = trace), method = "Nelder-Mead" , hessian = hess)
  }
  if (method == "SANN")  {
    w<-optim(init.par, fn = logL.joint.multi.R, C = C, y = y, m = m, n = n, anc.values = anc.values, yy = yy,
             control = list(fnscale = -1, maxit=10000, trace = trace), method = "SANN" , hessian = hess, lower = lower.limit)
  }

    if (w$convergence == 1) converge<-"Model did not converge"
    if (w$convergence == 0) converge<-"Model converged successfully"
    
  }
    # number of parameters
    K <- length(init.par) + (states-1) #parameters in the R matrices + ancestral values for each trait + number of shift points


  if (is.numeric(iterations) ==TRUE) {
  iter<-iterations
  }

  if (hess) {
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
    SE.R1<-matrix(0, nrow=m, ncol=m)
    SE.R2<-matrix(0, nrow=m, ncol=m)
    diag(SE.R1)<-w$se[1:m]
    diag(SE.R2)<-w$se[(m+1):(m+m)]

    locations.SE.R1<-which(SE.R1 == 0, arr.ind = T)
    location.upper.tri.R<-which(locations.SE.R1[,1] < locations.SE.R1[,2])

    upper.first<-w$se[(m+m+1):(m+m+length(location.upper.tri.R))]
    upper.second<-w$se[(m+m+length(location.upper.tri.R)+1):(m+m+length(location.upper.tri.R)+length(location.upper.tri.R))]

    for (i in 1:m){
      SE.R1[locations.SE.R1[,1][location.upper.tri.R[i]],locations.SE.R1[,2][location.upper.tri.R[i]]]<-upper.first[i]
      SE.R2[locations.SE.R1[,1][location.upper.tri.R[i]],locations.SE.R1[,2][location.upper.tri.R[i]]]<-upper.second[i]
    }

    SE.R.1<-t(SE.R1)%*%SE.R1
    SE.R.2<-t(SE.R2)%*%SE.R2
    SE.R <- list(SE.R.1 = SE.R.1,  SE.R.2 = SE.R.2)
    SE.anc <- w$se[((length(init.trait.var)*2 + length(init.cov.traits)*2)+1) : length(init.par)]
  }
  if (hess == FALSE) {
    SE.R<-NA
    SE.anc <-NA
  }

  chole.1<-matrix(0, nrow=m, ncol=m)
  chole.2<-matrix(0, nrow=m, ncol=m)
  diag(chole.1)<-w$par[1:m]
  diag(chole.2)<-w$par[(m+1):(m+m)]

  locations.R1<-which(chole.1 == 0, arr.ind = T)
  location.upper.tri.R<-which(locations.R1[,1] < locations.R1[,2])

  upper.first<-w$par[(m+m+1):(m+m+length(location.upper.tri.R))]
  upper.second<-w$par[(m+m+length(location.upper.tri.R)+1):(m+m+length(location.upper.tri.R)+length(location.upper.tri.R))]

  for (i in 1:m){
    chole.1[locations.R1[,1][location.upper.tri.R[i]],locations.R1[,2][location.upper.tri.R[i]]]<-upper.first[i]
    chole.2[locations.R1[,1][location.upper.tri.R[i]],locations.R1[,2][location.upper.tri.R[i]]]<-upper.second[i]
  }

  R.1<-t(chole.1)%*%chole.1
  R.2<-t(chole.2)%*%chole.2

  R <- list(R.1 = R.1,  R.2= R.2)

  if (is.numeric(iterations) == FALSE) {
    iter<-NA
  }

  ancestral.values<-w$par[((length(init.trait.var)*2 + length(init.cov.traits)*2)+1) : length(init.par)]

  wc<-as.evoTS.multi.BW.fit(converge, modelName = "Multivariate Random walk with two R matrices (with non-zero off-diagonal elements)", logL = w$value, ancestral.values = ancestral.values, SE.anc = SE.anc, R = R, SE.R = SE.R,
                              method = "Joint", K = K, n = length(yy$xx[,1]), iter=iter)
  return(wc)
}
