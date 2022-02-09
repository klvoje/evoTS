#'@title Class for fit to evolutionary sequence (time-series) models
#'
#'@description A function that combines useful information summarizing model fit.
#'
#'@param logL log-likelihood of model
#'
#'@param ancestral.values maximum-likelihood estimates of the ancestral trait values
#'
#'@param optima maximum-likelihood estimates of the optima
#'
#'@param A maximum-likelihood estimates of the parameters in the A matrix
#'
#'@param half.life the calculated half-life of the evolutionary process
#'
#'@param R maximum-likelihood estimates of the parameters in the R matrix
#'
#'@param method the parameterization used: Joint
#'
#'@param K number of parameters in the model
#'
#'@param n sample size
#'
#'@param iter the number of times the optimization method is run from different starting points. Default is NULL, meaning the optimization is run once.
#'
#'@details This function is used by the model-fitting routines for the multivariate Ornstein-Uhlenbeck models to create standardized output
#'
#'@note This function is not likely to be called directly by the user.
#'
#'@author Kjetil Lysne Voje
#'


as.evoTS.multi.OU.fit<-function (logL, ancestral.values, optima, A, half.life, R, method, K, n, iter)
{
  ic <- paleoTS::IC(logL = logL, K = K, n = n, method = "AICc")
  y <- list(logL = logL, AICc = ic, ancestral.values = ancestral.values, optima = optima, A = A, half.life = half.life, R = R,
            method = method, K = K, n = n, iter = iter)
  class(y) <- "paleoTSfit"
  return(y)
}

