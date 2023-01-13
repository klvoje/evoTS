#' @title Calculate the log-likelihood surface for a part of parameter space
#'
#' @description Function to calculate the log-likelihood surface for a part of parameter space for a General Random Walk.
#'
#' @param y an univariate paleoTS object.
#'
#'@param vstep.vec vector containing the parameter values of the variance parameter to be evaluated
#'
#'@param mstep.vec vector containing the parameter values of the mstep parameter to be evaluated
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
#' @importFrom plotly layout
#'
#'@examples
#'
#'## Generate a paleoTS objects
#'x <- paleoTS::sim.GRW(30)
#'
#'## Fit the the model to the data.
#'x1<-paleoTS::opt.joint.GRW(x)
#'
#'\donttest{
#'## Create log-likelihood surface (the example may take > 5 seconds to run)
#'loglik.surface.GRW(x, mstep.vec= seq(0,0.3,0.01), vstep.vec = seq(0,0.3,0.01))
#'}

loglik.surface.GRW<-function(y, mstep.vec, vstep.vec, pool = TRUE){

  if (min(vstep.vec)<0) stop("the vstep parameter cannot take negative values. The smallest value is 0")

  if (pool)
  y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  y$tt <- y$tt - min(y$tt)
  n <- length(y$mm)
  anc<-paleoTS::opt.joint.GRW(y)$parameters[1]
  loglik<-matrix(NA, ncol=length(mstep.vec), nrow=length(vstep.vec))
  mstep.limit<-rep(NA, length(mstep.vec))
  vstep.limit<-rep(NA, length(vstep.vec))

  for (i in 1:length(mstep.vec)){
    for (j in 1:length(vstep.vec)){
      mstep <- mstep.vec[i]
      vstep <- vstep.vec[j]
      VV <- vstep * outer(y$tt, y$tt, FUN = pmin)
      diag(VV) <- diag(VV) + y$vv/y$nn
      M <- rep(anc, n) + mstep * y$tt
      loglik[j,i] <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VV, log = TRUE)
    }
  }

  loglik<-loglik-max(loglik)
  loglik<-loglik+2
  loglik[loglik<0] <- 0

  for (i in 1:length(mstep.vec)){
    if (all(loglik[,i]==0)) mstep.limit[i]<-NA else mstep.limit[i]<-mstep.vec[i]
  }

  for (i in 1:length(vstep.vec)){
    if (all(loglik[i,]==0)) vstep.limit[i]<-NA else vstep.limit[i]<-vstep.vec[i]
  }

  out<-matrix(c(min(na.exclude(mstep.limit)), max(na.exclude(mstep.limit)), min(na.exclude(vstep.limit)), max(na.exclude(vstep.limit))), ncol=2, byrow = TRUE)
  colnames(out)<-c("lower", "upper")
  rownames(out)<-c("mstep", "vstep")
  print(out)

  plot_ly() %>% add_surface(x = ~mstep.vec, y = ~vstep.vec, z = ~loglik, type = 'mesh3d')%>%
  layout(
    scene = list(
      xaxis = list(title = "mstep"),
      yaxis = list(title = "vstep"),
      zaxis = list(title = "log-likelihood")
    )
  )


}

