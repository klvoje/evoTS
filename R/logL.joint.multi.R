#' @title Log-likelihoods for evolutionary models
#'
#' @description Returns log-likelihood for a multivariate Unbiased Random Walk model fitted to separate segments of a multivariate time series.
#'
#' @param init.par initial (starting) parameters values
#'
#' @param C distance matrix
#'
#' @param y vector containing all trait values from all traits
#'
#' @param m number of traits
#'
#' @param n number of populations
#'
#' @param anc.values initial values for the ancestral trait values
#'
#' @param yy a multivariate evoTS object
#'
#' @details In general, users will not be access these functions directly, but instead use the optimization functions, which use these functions to find the best-supported parameter values.
#'
#'@return The log-likelihood of the parameter estimates, given the data.
#'
#'@author Kjetil Lysne Voje


logL.joint.multi.R<-function (init.par, C , y , m , n, anc.values, yy)
{
  m<-length(anc.values)
  var.covar.1<-diag(c(rep(0,m)))
  var.covar.2<-diag(c(rep(0,m)))
  diag(var.covar.1)<-c(init.par[1:m])
  diag(var.covar.2)<-c(init.par[(m+1):(m+m)])
  locations.R<-which(var.covar.1 == 0, arr.ind = T)
  location.upper.tri.R<-which(locations.R[,1] < locations.R[,2])

  upper.first<-init.par[(m+m+1):(m+m+length(location.upper.tri.R))]
  upper.second<-init.par[(m+m+length(location.upper.tri.R)+1):(m+m+length(location.upper.tri.R)+length(location.upper.tri.R))]

  for (i in 1:m){
    var.covar.1[locations.R[,1][location.upper.tri.R[i]],locations.R[,2][location.upper.tri.R[i]]]<-upper.first[i]
    var.covar.2[locations.R[,1][location.upper.tri.R[i]],locations.R[,2][location.upper.tri.R[i]]]<-upper.second[i]
  }

  cholesky <- array( c( var.covar.1 , var.covar.2 ) , dim = c( m , m , 2 ) )

  M.init<-init.par[(m+m+length(location.upper.tri.R)+length(location.upper.tri.R)+1):(m+m+length(location.upper.tri.R)+length(location.upper.tri.R)+m)]
  M_temp<-matrix(data=NA, nrow=m, ncol=n)
  for (i in 1:m){
    M_temp[i,] <- rep(M.init[i], n)
  }
  M<-c(t(M_temp)) # vectorize M


    V <- matrix(0, nrow=length(M), ncol=length(M)) # making a variance-covariance matrix with dimensionality of n*m * n*m
    # computing V as the kronecker product of the Cholesky decomposed R matrices multiplied with distance matrices (C)
    tmp.1<-V + kronecker(t(cholesky[, , 1]) %*% cholesky[, , 1], C[[1]])
    tmp.2<-V + kronecker(t(cholesky[, , 2]) %*% cholesky[, , 2], C[[2]])
    VV<-tmp.1+tmp.2

    sample.var_temp<-matrix(data=NA, nrow=m, ncol=n)
    for (i in 1:m){
      sample.var_temp[i,] <- yy$vv[,i]/yy$nn[,i]
    }
    sample.var<-c(t(sample.var_temp))
    diag(VV) <- diag(VV) + sample.var

    y<-as.vector(y)

    S <- mvtnorm::dmvnorm(y, mean = M, sigma = VV, log = TRUE)
}
