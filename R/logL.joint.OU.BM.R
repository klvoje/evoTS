#' @title Log-likelihoods for evolutionary models
#'
#' @description Returns log-likelihood for an Ornstein-Uhlenbeck model where the optimum evolves as a Unbiased Random Walk. The movement of the optimum is not parameterized based on separate data.
#'
#' @param p parameters of the model to be optimized
#'
#' @param y a paleoTS object
#'
#' @param opt.anc logical, indicating if the ancestral trait value is at the optimum (TRUE) or displaced from the optimum (FALSE)
#'
#' @details In general, users will not be access these functions directly, but instead use the optimization functions, which use these functions to find the best-supported parameter values.
#'
#'@return The log-likelihood of the parameter estimates, given the data.
#'
#'@author Kjetil Lysne Voje
#'
#'@references Hansen, T. F., Pienaar, J. & Orzack, S. H. A Comparative Method for Studying Adaptation to a Randomly Evolving Environment. \emph{Evolution} 62, 1965–1977 (2008).

logL.joint.OU.BM<-function (p, y, opt.anc)
{

  if (opt.anc == TRUE){
    anc <- p[1]
    vs <- p[2]
    aa <- p[3]
    vo<-p[4]
    theta.0<-y$mm[1]
  }
  if (opt.anc == FALSE){
  anc <- p[1]
  vs <- p[2]
  theta.0 <- p[3]
  aa <- p[4]
  vo<-p[5]
  }

  n <- length(y$mm)
  ff <- function(a, b) abs(a - b)
  tij<-outer(y$tt, y$tt, FUN = ff)
  ta<-outer(y$tt, y$tt, pmin)
  #VCOV<- ((vo+vs)/(2 * aa)) * (1 - exp(-2 * aa * ta)) * exp(-aa *tij) + (vo*ta*(1-((1+exp(-aa*tij))*(1-exp(-aa*ta))/(aa*ta))))
  VCOV<- ((vo+vs)/(2 * aa)) * (1 - exp(-2 * aa * ta)) * exp(-aa *tij) + vo*ta*(1-(1+exp(-aa*tij))*(1-exp(-aa*ta))/(aa*ta))
  
  VCOV[1,]<-0;
  VCOV[,1]<-0;
  #VCOV<-((vo+vs)/(2 * aa)) * (1 - exp(-2 * aa * y$tt))   + (vo*y$tt*(1 - ((2*(1-exp(-aa*y$tt)))/(aa*y$tt))))
  diag(VCOV) <- diag(VCOV) + y$vv/y$nn
  theta<-rep(theta.0, n)
  M <- ou.M(anc, theta, aa, y$tt)
  S <- mvtnorm::dmvnorm(t(y$mm), mean = M, sigma = VCOV, log = TRUE)
  return(S)

}
