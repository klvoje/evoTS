#' @title Fit three OU models to three segments of an evolutionary sequence (time-series).
#'
#' @description Function to find maximum likelihood solutions to three OU models to three different segments of an evolutionary sequence. Most users will interact with this function through the function fit.OU.OU.OU.
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
#'@references Hunt, G., Bell, M. A. & Travis, M. P. Evolution toward a new adaptive optimum: Phenotypic evolution in a fossil stickleback lineage. \emph{Evolution} 62:700–710.
#'
#'@export

opt.joint.OU.OU.OU<-function (y, gg, pool = TRUE, meth = "L-BFGS-B", hess = FALSE) 
{
  if (pool) {
    y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  }
  
  cl = list(fnscale = -1)
  
  small <- 1e-08
  p0OU_1<- paleoTS::mle.GRW(paleoTS::sub.paleoTS(y, ok = gg == 1))
  p0OU_2<- paleoTS::mle.GRW(paleoTS::sub.paleoTS(y, ok = gg == 2))
  p0OU_3<- paleoTS::mle.GRW(paleoTS::sub.paleoTS(y, ok = gg == 3))
  
  
  K <- 12
  
  p0anc<-y$mm[1]
  names(p0anc)<-"anc"
  
  #prepare initial guesses of OU parameters 
  halft_1 <- (paleoTS::sub.paleoTS(y, ok = gg == 1)$tt[length(paleoTS::sub.paleoTS(y, ok = gg == 1)$tt)] - paleoTS::sub.paleoTS(y, ok = gg == 1)$tt[1])/4
  p0_1 <- c(p0OU_1["vstep"]/10, paleoTS::sub.paleoTS(y, ok = gg == 1)$mm[length(paleoTS::sub.paleoTS(y, ok = gg == 1)$mm)], log(2)/halft_1)
  names(p0_1) <- c("vstep_1", "theta_1", "alpha_1")
  
  halft_2 <- (paleoTS::sub.paleoTS(y, ok = gg == 2)$tt[length(paleoTS::sub.paleoTS(y, ok = gg == 2)$tt)] - paleoTS::sub.paleoTS(y, ok = gg == 2)$tt[1])/4
  p0_2 <- c(p0OU_2["vstep"]/10, paleoTS::sub.paleoTS(y, ok = gg == 2)$mm[length(paleoTS::sub.paleoTS(y, ok = gg == 2)$mm)], log(2)/halft_2)
  names(p0_2) <- c("vstep_2", "theta_2", "alpha_2")
  
  halft_3 <- (paleoTS::sub.paleoTS(y, ok = gg == 3)$tt[length(paleoTS::sub.paleoTS(y, ok = gg == 3)$tt)] - paleoTS::sub.paleoTS(y, ok = gg == 3)$tt[1])/4
  p0_3 <- c(p0OU_3["vstep"]/10, paleoTS::sub.paleoTS(y, ok = gg == 3)$mm[length(paleoTS::sub.paleoTS(y, ok = gg == 3)$mm)], log(2)/halft_2)
  names(p0_3) <- c("vstep_3", "theta_3", "alpha_3")
  
  #make vector out of estimated model parameters
  p0 <- c(p0anc, p0_1, p0_2, p0_3)
  
  #define lower bounds for the L-BFGS-B method
  # three parameters in the unbiased ramdom walk and 4 parameters for the biased random walk
  ll<- c(NA, small, NA, small, small, NA, small, small, NA, small)
  
  # If user wants to use a differen method than L-BFGS-B, then the lower boun is defined as -Inf
  if (meth != "L-BFGS-B") 
    ll <- -Inf
  
  #Here the multivariate parameter estimation routine start:
  
 w <- try(optim(p0, fn = logL.joint.OU.OU.OU, gg = gg, 
               method = meth, lower = ll, control = cl, hessian = hess, 
                 y = y), silent = TRUE)
  if (class(w) == "try-error") {
    cl <- list(fnscale = -1, parscale = c(1, 10, 100))
    w <- try(optim(p0, fn = logL.joint.OU.OU.OU, gg = gg, 
                   method = meth, lower = ll, control = cl, hessian = hess, 
                   y = y), silent = TRUE)
  }
  if (class(w) == "try-error") {
    wc <- paleoTS::as.paleoTSfit(logL = NA, parameters = NA, modelName = paste("OU", 
                                                                      "OU", "OU", sep = "-"), method = "Joint", K = K, n = length(y$mm), 
                        se = NULL)
    return(wc)
  }
  if (hess) 
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- paleoTS::as.paleoTSfit(logL = w$value, parameters = w$par, modelName = paste("OU", 
                                                                            "OU", "OU", sep = "-"), method = "Joint", K = K, n = length(y$mm), 
                      se = w$se)
  return(wc)
}
