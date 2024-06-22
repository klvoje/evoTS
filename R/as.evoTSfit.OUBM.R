#'@title Class for fit to evolutionary sequence (time-series) models
#'
#'@description A function that combines useful information summarizing model fit.
#'
#'@param logL log-likelihood of model
#'
#'@param parameters maximum-likelihood estimates of the parameters
#'
#'@param modelName description of the model
#'
#'@param method the parameterization used: Joint
#'
#'@param K number of parameters in the model
#'
#'@param n sample size
#'
#'@param se standard errors of parameter estimates
#'
#'@details This function is used by the model-fitting routines for the univariate Ornstein-Uhlenbeck model where the optimum evolves as an Unbiased Random Walk to create standardized output
#'
#'@note This function is not likely to be called directly by the user.
#'
#'@author Kjetil Lysne Voje
#'

as.evoTSfit.OUBM<-function (logL, parameters, modelName, method, K, n, se)
{
  ic <- paleoTS::IC(logL = logL, K = K, n = n, method = "AICc")
  y <- list(logL = logL, AICc = ic, parameters = parameters,
            modelName = modelName, method = method, K = K, n = n, se = se)
  class(y) <- "paleoTSfit"
  return(y)
}
