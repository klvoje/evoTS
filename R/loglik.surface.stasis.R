#' @title Calculate the log-likelihood surface for a part of parameter space
#'
#' @description Function to calculate the log-likelihood surface for a part of parameter space for the Stasis model.
#'
#' @param y an univariate paleoTS object.
#'
#'@param theta.vec vector containing the parameter values of the theta parameter to be evaluated
#'
#'@param omega.vec vector containing the parameter values of the omega parameter to be evaluated
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
#'x <- paleoTS::sim.Stasis(30)
#'
#'## Fit the the model to the data.
#'x1<-paleoTS::opt.joint.Stasis(x)
#'
#'\donttest{
#'## Create log-likelihood surface (the example may take > 5 seconds to run)
#'loglik.surface.stasis(x, theta.vec= seq(-0.15,0.1,0.001), omega.vec = seq(0,0.1,0.001))
#'}

loglik.surface.stasis<-function(y, theta.vec, omega.vec, pool = TRUE){

  x<-NULL

  if (min(omega.vec)<0) stop("the omega parameter cannot take negative values. The smallest value is 0")

  if (pool)
  y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  y$tt <- y$tt - min(y$tt)
  n <- length(y$mm)
  loglik<-matrix(NA, ncol=length(theta.vec), nrow=length(omega.vec))
  theta.limit<-rep(NA, length(theta.vec))
  omega.limit<-rep(NA, length(omega.vec))

  for (i in 1:length(theta.vec)){
    for (j in 1:length(omega.vec)){
      theta <- theta.vec[i]
      omega <- omega.vec[j]
      VV <- diag(omega + y$vv/y$nn)
      M <- rep(theta, n)
      loglik[j,i] <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VV, log = TRUE)
    }
  }

  loglik<-loglik-max(loglik)
  loglik<-loglik+2
  loglik[loglik<0] <- 0

  for (i in 1:length(theta.vec)){
    if (all(loglik[,i]==0)) theta.limit[i]<-NA else theta.limit[i]<-theta.vec[i]
  }

  for (i in 1:length(omega.vec)){
    if (all(loglik[i,]==0)) omega.limit[i]<-NA else omega.limit[i]<-omega.vec[i]
  }

  out<-matrix(c(min(na.exclude(theta.limit)), max(na.exclude(theta.limit)), min(na.exclude(omega.limit)), max(na.exclude(omega.limit))), ncol=2, byrow = TRUE)
  colnames(out)<-c("lower", "upper")
  rownames(out)<-c("theta", "omega")
  print(out)

  plot_ly() %>% add_surface(x = ~theta.vec, y = ~omega.vec, z = ~loglik, type = 'mesh3d')%>%
  layout(
    scene = list(
      xaxis = list(title = "theta"),
      yaxis = list(title = "omega"),
      zaxis = list(title = "log-likelihood")
    )
  )


}

