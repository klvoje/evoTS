#' @title Log-likelihoods for evolutionary models
#'
#' @description Returns log-likelihood for an Unbiased Random Walk with an accelerating or decelerating rate of change through time.
#'
#' @param p parameters of the model to be optimized
#'
#' @param y a paleoTS object
#'
#' @details In general, users will not be access these functions directly, but instead use the optimization functions, which use these functions to find the best-supported parameter values.
#'
#'@return The log-likelihood of the parameter estimates, given the data.
#'
#'@author Kjetil Lysne Voje


logL.joint.accel_decel<-function (p, y)
{
  anc <- p[1]
  vs <- p[2]
  r <- p[3]
  n <- length(y$mm)
  VV <- vs * outer(exp(r*y$tt), exp(r*y$tt), FUN = pmin)
  diag(VV) <- diag(VV) + y$vv/y$nn
  M <- rep(anc, n)
  S <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VV, log = TRUE)
  return(S)
}
