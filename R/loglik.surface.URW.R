#' @title Calculate the log-likelihood surface for a part of parameter space
#'
#' @description Function to calculate the log-likelihood surface for a part of parameter space for a Unbiased Random Walk.
#'
#' @param y an univariate paleoTS object.
#'
#'@param vstep.vec vector containing the parameter values of the variance parameter to be evaluated
#'
#'@param pool indicating whether to pool variances across samples
#'
#'@return the function returns the range of parameter values that are within two log-likelihood units from the best (maximum) parameter estimate and a log-likelihood surface.
#'
#'@note How fine-scaled the estimated log-likelihood surface is depends on the step size between the values in the input-vectors. The step-size therefore determines how accurate the representation of the support surface is, including the returned upper and lower estimates printed in the console. The range of the input vectors needs to be increased if the confidence interval includes the boundary of the input vector. Note also that it might be wise to include the maximum likelihood estimates as part of the input vectors. The computed support surface is conditional on the best estimates of the other model parameters that are not part of the support surface (e.g. the estimated ancestral trait value).
#'
#'@author Kjetil Lysne Voje
#'
#'@export
#'
#'@examples
#'
#'## Generate a paleoTS objects
#'x <- paleoTS::sim.GRW(30)
#'
#'## Fit a the model to the data by defining shift points.
#'x1<-paleoTS::opt.joint.URW(x)
#'
#'## Create log-likelihood surface (the example may take > 5 seconds to run)
#'loglik.surface.URW(x, vstep.vec = seq(0,0.5,0.001))
#'

loglik.surface.URW<-function(y, vstep.vec, pool = TRUE){

  x<-NULL
  if (min(vstep.vec)<0) stop("the vstep parameter cannot take negative values. The smallest value is 0")

  y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  y$tt <- y$tt - min(y$tt)
  n <- length(y$mm)
  anc<-paleoTS::opt.joint.URW(y)$parameters[1]
  loglik<-rep(NA, length(vstep.vec))
  vstep.limit<-rep(NA, length(vstep.vec))

    for (i in 1:length(vstep.vec)){
      vstep <- vstep.vec[i]
      VV <- vstep * outer(y$tt, y$tt, FUN = pmin)
      diag(VV) <- diag(VV) + y$vv/y$nn
      M <- rep(anc, n)
      loglik[i] <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VV, log = TRUE)
    }

  loglik<-loglik-max(loglik)
  loglik<-loglik+2
  loglik[loglik<0] <- 0

  for (i in 1:length(vstep.vec)){
    if (all(loglik[i]==0)) vstep.limit[i]<-NA else vstep.limit[i]<-vstep.vec[i]
  }

  out<-matrix(c(min(na.exclude(vstep.limit)), max(na.exclude(vstep.limit))), ncol=2, byrow = TRUE)
  colnames(out)<-c("lower", "upper")
  rownames(out)<-"vstep"
  print(out)
  plot(loglik~vstep.vec, type="l", col="black", lwd=3, xlab="vstep", ylab="log-likelihood", cex.lab=1.2)

}


