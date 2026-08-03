#' @title Fit multivariate Unbiased Random Walk model.
#'
#' @description Function to find maximum likelihood solution to a multivariate Unbiased Random Walk model.
#'
#' @param yy a multivariate evoTS object.
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
#' @details The function searches - using an optimization routine - for the maximum-likelihood solution for a multivariate Unbiased Random Walk model.
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
#'## Fit a multivariate Unbiased Random Walk model.
#'opt.single.R(x)

opt.single.R<-function (yy, method="L-BFGS-B", hess = FALSE, pool = TRUE, trace=FALSE, iterations=NULL, iter.sd=NULL)
{

  n <- nrow(yy$xx) #number of samples/populations
  m <- ncol(yy$xx) # number of traits
  
  if (pool==TRUE) { 
    for (i in 1:m){
      
      tmp<-paleoTS::as.paleoTS(yy$xx[,i], yy$vv[,i], yy$nn[,i], yy$tt[,i])
      tmp<- paleoTS::pool.var(tmp, ret.paleoTS = TRUE)
      yy$vv[,i]<-tmp$vv
    }
  }

  C<- outer(yy$tt[,1], yy$tt[,1], FUN = pmin) #Create distance matrix
### v.1.4 ###
# sampling-error vector computed once here instead of being rebuilt from
# yy inside the log-likelihood function on every optimizer call.
  se_vec <- as.vector(t(yy$vv / yy$nn))  # pre-compute sampling error vector

  X <- yy$xx # Character matrix with dimensions n * m
  y <- as.matrix(as.vector(X)) # Vectorized version of X


  # Define initial parameter values for the optimization routine
  init.trait.var<-apply(yy$xx,2,var)
### v.1.4 ###
# V.1.3 used cov(as.matrix(yy$xx)) (covariance of raw trait values)
# as the initial trait-covariance guess. Now uses the covariance of the
# differenced series, a more appropriate starting guess for a rate parameter.
  temp.matrix<-cov(apply(yy$xx, 2, diff))
  init.cov.traits<-unique(temp.matrix[row(temp.matrix)!=col(temp.matrix)])
  anc.values<-yy$xx[1,]

  init.par<-c(init.trait.var, init.cov.traits, anc.values)
  lower.limit<-c(rep(0,length(init.trait.var)), rep(NA,length(init.cov.traits)), rep(NA, length(anc.values)))
### v.1.4 ###
# per-trait data SD, used below to scale additive perturbations of
# ancestral values during random restarts.
  data.sd <- apply(yy$xx, 2, stats::sd)

  ### Start iterations from different starting values
  if (is.numeric(iterations)) {
    if(is.numeric(iter.sd) == FALSE) iter.sd <-1
    #if(is.numeric(max.attemps) == FALSE) max.attemps <-100000
    log.lik.tmp<-rep(NA, 1000000)
    www<-list()

    for (k in 1:1000000){
      tryCatch({
### v.1.4 ###
# Verion 1.3 perturbed the whole parameter vector identically:
# `init.par_temp<-init.par; init.par<-rnorm(length(init.par_temp), init.par_temp, iter.sd)`.
# Now each parameter type is perturbed in a scale-appropriate way: multiplicative
# log-normal jitter for the (positive) trait-variance parameters, small additive
# jitter for the covariance terms, and data-SD-scaled additive jitter for the
# ancestral values. Result is also clamped to the lower bound (NA bounds treated
# as -Inf) so restarts can't violate constraints.
      init.trait.var  <- init.trait.var * exp(rnorm(length(init.trait.var), 0, iter.sd))
      init.cov.traits <- rnorm(length(init.cov.traits), init.cov.traits, iter.sd * 0.1)
      anc.values      <- anc.values + rnorm(length(anc.values), 0, data.sd * iter.sd)
      init.par        <- c(init.trait.var, init.cov.traits, anc.values)
      init.par        <- pmax(replace(lower.limit, is.na(lower.limit), -Inf), init.par)

      if (method == "L-BFGS-B")  {
### v.1.4 ###
# `yy = yy` replaced by `se_vec = se_vec` (see pre-computation above); the
# `lower` bound now replaces NA entries with -Inf before passing to optim(),
# since optim() itself does not accept NA bounds.
        www[[k]]<-optim(init.par, fn = logL.joint.single.R, C = C, y = y, m = m, n = n, anc.values = anc.values, se_vec = se_vec,
                        control = list(fnscale = -1, maxit=10000, trace = trace), method = "L-BFGS-B", hessian = hess, lower = replace(lower.limit, is.na(lower.limit), -Inf))
      }

      if (method == "Nelder-Mead")  {
### v.1.4 ###
# `yy = yy` replaced by `se_vec = se_vec` (see pre-computation above).
        www[[k]]<-optim(init.par, fn = logL.joint.single.R, C = C, y = y, m = m, n = n, anc.values = anc.values, se_vec = se_vec,
                        control = list(fnscale = -1, maxit=10000, trace = trace), method = "Nelder-Mead" , hessian = hess)
      }
      if (method == "SANN")  {
### v.1.4 ###
# `yy = yy` replaced by `se_vec = se_vec` (see pre-computation above).
        www[[k]]<-optim(init.par, fn = logL.joint.single.R, C = C, y = y, m = m, n = n, anc.values = anc.values, se_vec = se_vec,
                        control = list(fnscale = -1, maxit=10000, trace = trace), method = "SANN" , hessian = hess, lower = lower.limit)
      }
      log.lik.tmp[k]<-www[[k]]$value
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      
      if (length(na.exclude(log.lik.tmp)) == iterations){
        break
      }
    }
    
    # Need to remove entries in www in case there are iterations where the initial parameter estimates did not work.
    www_tmp<-list()
    for (i in 1:k){
      if (is.character(www[[i]][1]) == FALSE) www_tmp[i]<-list(www[[i]])
    }
    
    www_reduced<-www_tmp[!sapply(www_tmp,is.null)]
    
    
    for (j in 1:length(www_reduced)){
      if(max(na.exclude(log.lik.tmp)) == www_reduced[[j]]$value) best.run<-www_reduced[[j]]
    }
    
  }


  ##### Start of non-iteration routine #####

  if (is.numeric(iterations) == FALSE) {


    if (method == "L-BFGS-B")  {
### v.1.4 ###
# `yy = yy` replaced by `se_vec = se_vec` (see pre-computation above); the
# `lower` bound now replaces NA entries with -Inf before passing to optim().
      w<-optim(init.par, fn = logL.joint.single.R, C = C, y = y, m = m, n = n, anc.values = anc.values, se_vec = se_vec,
               control = list(fnscale = -1, maxit=10000, trace = trace), method = "L-BFGS-B", hessian = hess, lower = replace(lower.limit, is.na(lower.limit), -Inf))
    }

    if (method == "Nelder-Mead")  {
### v.1.4 ###
# `yy = yy` replaced by `se_vec = se_vec` (see pre-computation above).
      w<-optim(init.par, fn = logL.joint.single.R, C = C, y = y, m = m, n = n, anc.values = anc.values, se_vec = se_vec,
               control = list(fnscale = -1, maxit=10000, trace = trace), method = "Nelder-Mead" , hessian = hess)
    }
    if (method == "SANN")  {
### v.1.4 ###
# `yy = yy` replaced by `se_vec = se_vec` (see pre-computation above).
      w<-optim(init.par, fn = logL.joint.single.R, C = C, y = y, m = m, n = n, anc.values = anc.values, se_vec = se_vec,
               control = list(fnscale = -1, maxit=10000, trace = trace), method = "SANN" , hessian = hess, lower = lower.limit)
    }
    
  }
  # number of parameters
  K <- length(init.par) #parameters in the R matrices + ancestral values for each trait

  if (is.numeric(iterations) ==TRUE) {
    iter<-iterations
    w<-best.run
  }
  
  if (w$convergence == 10) converge<-"The search algorithm stopped as it did not make progress towards the optimal solution"
  if (w$convergence == 0) converge<-"Model converged successfully"
  if (w$convergence == 1) converge<-"The maximum number of iterations was reached and the search algorithm exited"
  if (w$convergence == 51) converge<-"The model did not converge due to en error in L-BFGS-B. Reported estimates are not the maximum likelihood"
  if (w$convergence == 52) converge<-"The model did not converge due to en error in L-BFGS-B. Reported estimates are not the maximum likelihood"
  

  if (hess) {
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
    SE.R1<-matrix(0, nrow=m, ncol=m)
    diag(SE.R1)<-w$se[1:m]

    locations.SE.R1<-which(SE.R1 == 0, arr.ind = T)
    location.upper.tri.R<-which(locations.SE.R1[,1] < locations.SE.R1[,2])

    upper.first<-w$se[(m+1):(m+length(location.upper.tri.R))]

    for (i in 1:m){
      SE.R1[locations.SE.R1[,1][location.upper.tri.R[i]],locations.SE.R1[,2][location.upper.tri.R[i]]]<-upper.first[i]
    }

    SE.R<-t(SE.R1)%*%SE.R1
    SE.anc<-tail(w$se,m)
  }

  if (hess == FALSE) {
    SE.R<-NA
    SE.anc<-NA
  }

  chole.1<-matrix(0, nrow=m, ncol=m)
  diag(chole.1)<-w$par[1:m]

  locations.R1<-which(chole.1 == 0, arr.ind = T)
  location.upper.tri.R<-which(locations.R1[,1] < locations.R1[,2])

  upper.first<-w$par[(m+1):(m+length(location.upper.tri.R))]

  for (i in 1:m){
    chole.1[locations.R1[,1][location.upper.tri.R[i]],locations.R1[,2][location.upper.tri.R[i]]]<-upper.first[i]
  }

  R<-t(chole.1)%*%chole.1

  if (is.numeric(iterations) == FALSE) {
    iter<-NA
  }

  ancestral.values<-w$par[(length(init.trait.var) + length(init.cov.traits) +1) : length(init.par)]

  wc<-as.evoTS.multi.URW.fit(converge, modelName = "Multivariate model: Random walk (R matrix with non-zero off-diagonal elements)", logL = w$value, ancestral.values = ancestral.values, SE.anc = SE.anc, R = R, SE.R = SE.R,
                                                          method = "Joint", K = K, n = length(yy$xx[,1]), iter=iter)

  return(wc)
}
