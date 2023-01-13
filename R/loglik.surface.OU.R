#' @title Calculate the log-likelihood surface for a part of parameter space
#'
#' @description Function to calculate the log-likelihood surface for a part of parameter space for a Ornstein-Uhlenbeck model.
#'
#' @param y an univariate paleoTS object.
#'
#'@param stat.var.vec vector containing the parameter values of the stationary variance to be evaluated
#'
#'@param h.vec vector containing the parameter values of the half life to be evaluated
#'
#'@param anc the ancestral state
#'
#'@param theta the optimum
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
#'x <- paleoTS::sim.OU(40)
#'
#'## Fit the model to the data.
#'x1<-paleoTS::opt.joint.OU(x)
#'
#'##calculate half-life from model output
#'log(2)/x1$parameters[4]
#'
#'##calculate stationary variance from model output
#'x1$parameters[2]/(2*x1$parameters[4])
#'
#'\donttest{
#'## Create log-likelihood surface (the example may take > 5 seconds to run)
#'loglik.surface.OU(x, stat.var.vec=seq(0.001,0.5,0.01), h.vec=seq(0.01,10, 0.1))
#'}


loglik.surface.OU<-function(y, stat.var.vec, h.vec, anc = NULL, theta = NULL, pool = TRUE){

  if (min(stat.var.vec)<0) stop("the stationary variance cannot take negative values. The smallest value is 0")
  if (min(h.vec)<0) stop("the half-life cannot take negative values. The smallest value is 0")

  if (pool)
  y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  y$tt <- y$tt - min(y$tt)
  n <- length(y$mm)
  if (is.numeric(anc) == TRUE) anc<-anc else anc <- paleoTS::opt.joint.OU(y)$parameters[1]
  if (is.numeric(theta) == TRUE) theta<-theta else theta <- theta<-paleoTS::opt.joint.OU(y)$parameters[3]


  loglik<-matrix(NA, ncol=length(stat.var.vec), nrow=length(h.vec))
  stat.var.limit<-rep(NA, length(stat.var.vec))
  h.limit<-rep(NA, length(h.vec))

  for (i in 1:length(stat.var.vec)){
    for (j in 1:length(h.vec)){
      aa <- log(2)/h.vec[j]
      vs <- stat.var.vec[i]*2*aa
      ff <- function(a, b) abs(a - b)
      VV <- outer(y$tt, y$tt, FUN = ff)
      VV <- exp(-aa * VV)
      VVd <- ou.V(vs, aa, y$tt)
      VV2 <- outer(VVd, VVd, pmin)
      VV <- VV * VV2
      diag(VV) <- VVd + y$vv/y$nn
      M <- ou.M(anc, theta, aa, y$tt)
      loglik[j,i] <- mvtnorm::dmvnorm(t(y$mm), mean = M, sigma = VV, log = TRUE)
    }
  }
  loglik<-loglik-max(loglik)
  loglik<-loglik+2
  loglik[loglik<0] <- 0

  for (i in 1:length(stat.var.vec)){
    if (all(loglik[,i]==0)) stat.var.limit[i]<-NA else stat.var.limit[i]<-stat.var.vec[i]
  }

  for (i in 1:length(h.vec)){
    if (all(loglik[i,]==0)) h.limit[i]<-NA else h.limit[i]<-h.vec[i]
  }

  out<-matrix(c(min(na.exclude(stat.var.limit)), max(na.exclude(stat.var.limit)), min(na.exclude(h.limit)), max(na.exclude(h.limit))), ncol=2, byrow = TRUE)
  colnames(out)<-c("lower", "upper")
  rownames(out)<-c("stationary variance","half-life")
  print(out)

  plot_ly() %>% add_surface(x = ~stat.var.vec, y = ~h.vec, z = ~loglik, type = 'mesh3d')%>%
  layout(
    scene = list(
      xaxis = list(title = "stationary variance"),
      yaxis = list(title = "half-life"),
      zaxis = list(title = "log-likelihood")
    )
  )


}

