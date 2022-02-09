#' @title Fit thee models (Stasis, General Random Walk (trend), Unbiased Random Walk) to three segments of an evolutionary sequence (time-series).
#'
#' @description Function to find maximum likelihood solutions to three models to three different segments of an evolutionary sequence. Most users will interact with this function through the function fit.Stasis.trend.Stasis.
#'
#' @param y an univariate evoTS object.
#' 
#' @param gg numeric vector indicating membership of each sample in segments.
#'
#' @param pool logical indicating whether to pool variances across samples
#'
#' @param meth optimization method, passed to function optim. Default is "L-BFGS-B".
#' 
#'@param hess logical, indicating whether to calculate standard errors from the Hessian matrix. 
#' 
#'@return 
#'\item{logL}{the log-likelihood of the optimal solution} 
#'\item{AICc}{AIC with a correction for small sample sizes}
#'\item{parameters}{parameter estimates}
#'\item{modelName}{abbreviated model name}
#'\item{method}{Joint consideration of all samples} 
#'\item{K}{number of parameters in the model} 
#'\item{n}{the number of observations/samples} 
#'\item{all.logl}{log-likelihoods for all tested partitions of the series into segments. Will return a single value if shift points have been given} 
#'\item{GG}{matrix of indices of initial samples of each tested segment configuration; each column of GG corresponds to the elements of all.logl} 
#'
#'@note The models have been implemented to be compatible with the joint parameterization routine in the package paleoTS. The optimization is therefore fit using the actual sample values, with the autocorrelation among samples accounted for in the log-likelihood function. The joint distribution of sample means is multivariate normal, with means and variance-covariances determined by evolutionary parameters and sampling errors. 
#'
#'@author Kjetil Lysne Voje 
#'
#'@references Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond. \emph{Paleobiology} 32:578–601. 
#'
#'@export

opt.joint.Stasis.trend.URW<-function (y, gg, pool = TRUE, meth = "L-BFGS-B", hess = FALSE) 
{
  if (pool) {
    y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  }
  
  cl = list(fnscale = -1)
  small <- 1e-08
    p0st_1 <- paleoTS::mle.Stasis(paleoTS::sub.paleoTS(y, ok = gg == 1))
    p0dt<- paleoTS::mle.GRW(paleoTS::sub.paleoTS(y, ok = gg == 2))
    p0rw<- paleoTS::mle.URW(paleoTS::sub.paleoTS(y, ok = gg == 3))
    
    K <- 7
  
  #Checking if te variance parameter in the RW and Stasis is below critical value, and if TRUE, multiply it with 100
  if (p0st_1["omega"] <= small) 
    p0st_1["omega"] <- 100 * small
  if (p0dt["vstep"] <= small) 
    p0dt["vstep"] <- 100 * small
  if (p0rw["vstep"] <= small) 
    p0rw["vstep"] <- 100 * small
  
    #make vector out of estimated model parameters
  p0 <- c(p0st_1, p0dt, p0rw)
  
  #define lower bounds for the L-BFGS-B method
  # three parameters in the unbiased ramdom walk and 4 parameters for the biased random walk
  ll<- c(NA, small, NA, small, small)
  
  # If user wants to use a differen method than L-BFGS-B, then the lower boun is defined as -Inf
  if (meth != "L-BFGS-B") 
    ll <- -Inf
  
    #Here the multivariate parameter estimation routine start:
    w <- try(optim(p0, fn = logL.joint.Stasis.trend.URW, gg = gg, 
                   method = meth, lower = ll, control = cl, hessian = hess, 
                   y = y), silent = TRUE)
    if (class(w) == "try-error") {
      cl <- list(fnscale = -1, parscale = c(1, 10, 100))
      w <- try(optim(p0, fn = logL.joint.Stasis.trend.URW, gg = gg, 
                     method = meth, lower = ll, control = cl, hessian = hess, 
                     y = y), silent = TRUE)
    }
  
  if (class(w) == "try-error") {
    wc <- paleoTS::as.paleoTSfit(logL = NA, parameters = NA, modelName = paste("Stasis", 
                                                                      "Trend","URW", sep = "-"), method = "Joint", K = K, n = length(y$mm), 
                        se = NULL)
    return(wc)
  }
  if (hess) 
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- paleoTS::as.paleoTSfit(logL = w$value, parameters = w$par, modelName = paste("Stasis", 
                                                                            "Trend","URW", sep = "-"), method = "Joint", K = K, n = length(y$mm), 
                      se = w$se)
  return(wc)
}