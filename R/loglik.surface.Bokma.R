#' @title Calculate the log-likelihood surface for a part of parameter space
#'
#' @description Function to calculate the log-likelihood surface for a part of parameter space for the Bokma model.
#'
#' @param y an univariate paleoTS object.
#'
#'@param vstep_1.vec vector containing the parameter values of the vstep parameter before the punctuation to be evaluated
#'
#'@param vstep_2.vec vector containing the parameter values of the vstep parameter after the punctuation to be evaluated
#'
#'@param shift the sample where the punctuation happen
#'
#'@param anc the assumed ancestral state
#'
#'@param new.opt the new optimum (ancestral state for the Unbiased Random Walk after the punctuation)
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
#'## Generate a paleoTS objects
#'x <- sim.Bokma()
#'
#'## Fit a the model to the data by defining shift points.
#'x1<-fit.Bokma(x, minb=10)
#'
#'v_1<-seq(0,1,0.01)
#'
#'v_2<-seq(0,1, 0.01)
#'
#'## Create log-likelihood surface (the example may take > 5 seconds to run)
#'## loglik.surface.Bokma(x, vstep_1.vec=v_1, vstep_2.vec=v_2, shift=11, anc=0.4280, new.opt=1.7025)

loglik.surface.Bokma<-function(y, vstep_1.vec, vstep_2.vec, shift, anc, new.opt, pool = TRUE){

  if (min(vstep_1.vec)<0) stop("the vstep parameter cannot take negative values. The smallest value is 0")
  if (min(vstep_2.vec)<0) stop("the vstep parameter cannot take negative values. The smallest value is 0")

  if (pool)
  y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  y$tt <- y$tt - min(y$tt)
  n <- length(y$mm)
  loglik<-matrix(NA, ncol=length(vstep_1.vec), nrow=length(vstep_2.vec))
  vstep_1.limit<-rep(NA, length(vstep_1.vec))
  vstep_2.limit<-rep(NA, length(vstep_2.vec))

  for (i in 1:length(vstep_1.vec)){
    for (j in 1:length(vstep_2.vec)){
      vstep_1 <- vstep_1.vec[i]
      vstep_2 <- vstep_2.vec[j]
      n <- length(y$mm)
      n.1 <- shift
      n.2 <- length(y$mm)-n.1
      #creates a vector that repeates the ancestral value (first sample mean in sequence)
      M.1 <- rep(anc, n.1)
      M.2 <- rep(new.opt, n.2)
      M.1 <- unname(M.1)
      M.2 <- unname(M.2)
      M<-c(M.1, M.2)

      gg <- shift2gg(shift, length(y$mm))
      #define a vector based on length of first and second RW in the time series
      rw_1.seg <- which(gg == 1)
      rw_2.seg <- which(gg == 2)

      #Extract the time vector for the random walks:
      tt.rw_1 <- y$tt[rw_1.seg]
      tt.rw_2 <- y$tt[rw_2.seg]
      tt.rw_2<-tt.rw_2-tt.rw_2[1]

      #Multiply variance parameter (vstep) from RW with outer product of time*time
      VVrw_1 <- vstep_1 * outer(tt.rw_1, tt.rw_1, FUN = pmin)
      VVrw_2 <- vstep_2 * outer(tt.rw_2, tt.rw_2, FUN = pmin)

      #Create empty n*n matrix
      VVtot <- array(0, dim = c(n, n))

      #Create final variance matrix by combining the variance matrices from stasis and RW
      VVtot[rw_1.seg, rw_1.seg] <- VVrw_1
      VVtot[rw_2.seg, rw_2.seg] <- VVrw_2

      #Add population variance to the diagonal of the variance matrix
      diag(VVtot) <- diag(VVtot) + y$vv/y$nn

      loglik[j,i] <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VVtot, log = TRUE)
    }
  }

  loglik<-loglik-max(loglik)
  loglik<-loglik+2
  loglik[loglik<0] <- 0

  for (i in 1:length(vstep_1.vec)){
    if (all(loglik[,i]==0)) vstep_1.limit[i]<-NA else vstep_1.limit[i]<-vstep_1.vec[i]
  }

  for (i in 1:length(vstep_2.vec)){
    if (all(loglik[i,]==0)) vstep_2.limit[i]<-NA else vstep_2.limit[i]<-vstep_2.vec[i]
  }

  out<-matrix(c(min(na.exclude(vstep_1.limit)), max(na.exclude(vstep_1.limit)), min(na.exclude(vstep_2.limit)), max(na.exclude(vstep_2.limit))), ncol=2, byrow = TRUE)
  colnames(out)<-c("lower", "upper")
  rownames(out)<-c("vstep_1", "vstep_2")
  print(out)

  plot_ly() %>% add_surface(x = ~vstep_1.vec, y = ~vstep_2.vec, z = ~loglik, type = 'mesh3d')%>%
  layout(
    scene = list(
      xaxis = list(title = "vstep_1"),
      yaxis = list(title = "vstep_2"),
      zaxis = list(title = "log-likelihood")
    )
  )


}

