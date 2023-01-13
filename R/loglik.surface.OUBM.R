#' @title Calculate the log-likelihood surface for a part of parameter space
#'
#' @description Function to calculate the log-likelihood surface for a part of parameter space for a Ornstein-Uhlenbeck model where the optimum changes according to an Unbiased Random Walk.
#'
#' @param y an univariate paleoTS object.
#'
#'@param stat.var.vec vector containing the parameter values of the stationary variance to be evaluated
#'
#'@param h.vec vector containing the parameter values of the half life to be evaluated
#'
#'@param anc the ancestral state
#'
#'@param theta.0 the optimum
#'
#'@param vo the variance (vstep) parameter for the optimum
#'
#'@param opt.anc logical, indicating whether the the ancestral trait state is at the optimum.
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
#'x <- sim.OUBM(40)
#'
#'## Fit the model.
#'x1<-opt.joint.OUBM(x)
#'
#'##calculate half-life from model output
#'log(2)/x1$parameters[3]
#'
#'##calculate stationary variance from model output
#'x1$parameters[2]/(2*x1$parameters[3])
#'
#'\donttest{
#'## Create log-likelihood surface (the example may take > 5 seconds to run)
#'loglik.surface.OUBM(x, stat.var.vec=seq(0,4,0.01), h.vec=seq(0.0,5, 0.1))
#'}


loglik.surface.OUBM<-function(y, stat.var.vec, h.vec, anc = NULL, theta.0 = NULL, vo=NULL, opt.anc=TRUE, pool = TRUE){

  x<-NULL

  if (min(stat.var.vec)<0) stop("the stationary variance cannot take negative values. The smallest value is 0")
  if (min(h.vec)<0) stop("the half-life cannot take negative values. The smallest value is 0")

  if (pool)
  y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  y$tt <- y$tt - min(y$tt)
  n <- length(y$mm)
  
  if(opt.anc==TRUE){
    anc <- opt.joint.OUBM(y)$parameters[1]
    theta.0 <- opt.joint.OUBM(y)$parameters[1]
    vo <- opt.joint.OUBM(y)$parameters[4]
  } 
  
  if(opt.anc==FALSE){
    anc <- opt.joint.OUBM(y)$parameters[1]
    theta.0 <- opt.joint.OUBM(y)$parameters[3]
    vo <- opt.joint.OUBM(y)$parameters[5]
  } 
  
  if (is.numeric(anc) == TRUE) anc<-anc 
  if (is.numeric(theta.0) == TRUE) theta.0<-theta.0 
  if (is.numeric(vo) == TRUE) vo<-vo 

  loglik<-matrix(NA, ncol=length(stat.var.vec), nrow=length(h.vec))
  stat.var.limit<-rep(NA, length(stat.var.vec))
  h.limit<-rep(NA, length(h.vec))

  for (i in 1:length(stat.var.vec)){
    for (j in 1:length(h.vec)){
      aa <- log(2)/h.vec[j]
      vs <- stat.var.vec[i]*2*aa
      ff <- function(a, b) abs(a - b)
      tij<-outer(y$tt, y$tt, FUN = ff)
      ta<-outer(y$tt, y$tt, pmin)
      VCOV<- ((vo+vs)/(2 * aa)) * (1 - exp(-2 * aa * ta)) * exp(-aa *tij) + (vo*ta*(1-((1+exp(-aa*tij))*(1-exp(-aa*ta))/(aa*ta))))
      VCOV[1,]<-0;
      VCOV[,1]<-0;
      diag(VCOV) <- diag(VCOV) + y$vv/y$nn
      theta<-rep(theta.0, n)
      M <- ou.M(anc, theta, aa, y$tt)
      loglik[j,i] <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VCOV, log = TRUE)

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

