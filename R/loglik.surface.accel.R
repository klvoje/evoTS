#' @title Calculate the log-likelihood surface for a part of parameter space
#'
#' @description Function to calculate the log-likelihood surface for a part of parameter space for an Unbiased Random Walk with an accelerated rate of evolution.
#'
#' @param y an univariate paleoTS object.
#'
#'@param vstep.vec vector containing the parameter values of the variance parameter to be evaluated
#'
#'@param r.vec vector containing the parameter values of the r parameter to be evaluated
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
#' @importFrom plotly plot_ly
#' @importFrom plotly add_surface
#' @importFrom plotly  %>%
#'
#'@examples
#'## Generate a paleoTS objects
#'x <- sim.accel.decel(50)
#'
#'## Fit the model to the data.
#'x1<-opt.joint.accel(x)
#'
#'\donttest{
#'## Create log-likelihood surface (the example may take > 5 seconds to run)
#'loglik.surface.accel(x, vstep.vec = seq(0,4,0.005), r.vec = seq(0.15,0.25,0.005))
#'}

loglik.surface.accel<-function(y, vstep.vec, r.vec, pool = TRUE){

  if (min(vstep.vec)<0) stop("the variance parameter cannot take negative values. The smallest value is 0")
  if (min(r.vec)<0) stop("the exponential decay paramneter cannot take negative values. The largest value is 0")

  if (pool)
  y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  y$tt <- y$tt - min(y$tt)
  n <- length(y$mm)
  anc<-opt.joint.accel(y)$parameters[1]
  loglik<-matrix(NA, ncol=length(vstep.vec), nrow=length(r.vec))
  vstep.limit<-rep(NA, length(vstep.vec))
  r.limit<-rep(NA, length(r.vec))

  for (i in 1:length(vstep.vec)){
    for (j in 1:length(r.vec)){
      vstep <- vstep.vec[i]
      r <- r.vec[j]
      VV <- vstep * outer(exp(r*y$tt), exp(r*y$tt), FUN = pmin)
      diag(VV) <- diag(VV) + y$vv/y$nn
      M <- rep(anc, n)
      loglik[j,i] <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VV, log = TRUE)
    }
  }

  loglik<-loglik-max(loglik)
  loglik<-loglik+2
  loglik[loglik<0] <- 0

  for (i in 1:length(vstep.vec)){
    if (all(loglik[,i]==0)) vstep.limit[i]<-NA else vstep.limit[i]<-vstep.vec[i]
  }

  for (i in 1:length(r.vec)){
    if (all(loglik[i,]==0)) r.limit[i]<-NA else r.limit[i]<-r.vec[i]
  }

  out<-matrix(c(min(na.exclude(vstep.limit)), max(na.exclude(vstep.limit)), min(na.exclude(r.limit)), max(na.exclude(r.limit))), ncol=2, byrow = TRUE)
  colnames(out)<-c("lower", "upper")
  rownames(out)<-c("vstep", "r")
  print(out)

  plot_ly() %>% add_surface(x = ~vstep.vec, y = ~r.vec, z = ~loglik, type = 'mesh3d')%>%
    layout(
    scene = list(
      xaxis = list(title = "vstep"),
      yaxis = list(title = "r"),
      zaxis = list(title = "log-likelihood")
    )
  )


}

