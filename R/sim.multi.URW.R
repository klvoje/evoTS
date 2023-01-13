#' @title Simulate multivariate evolutionary sequence data that evolve according to an Unbiased Random Walk
#'
#' @description Function to simulate multivariate evolutionary sequence data that evolve according to an Unbiased Random Walk
#'
#' @param ns number of samples in time-series
#'
#' @param anc the ancestral trait values
#'
#' @param R the drift matrix
#'
#' @param vp within-population trait variance
#'
#' @param nn 	vector of the number of individuals in each sample (identical sample sizes for all time-series is assumed)
#'
#' @param tt 	vector of sample ages, increases from oldest to youngest
#'
#'@return A multivariate evolutionary sequence (time-series) data set.
#'
#'@author Kjetil Lysne Voje
#'
#'@export
#'
#'@examples
#'## Create a multivariate dataset
#'data_set<-sim.multi.URW(40, R = matrix(c(0.2,0.1,0.1,0.3), nrow=2, byrow = TRUE))
#'
#'## plot the data
#'plotevoTS.multivariate(data_set)
#'

sim.multi.URW<-function(ns = 30, anc = c(0,0), R = matrix(c(0.5,0,0,0.5), nrow=2, byrow = TRUE),
                       vp = 0.1, nn = rep(30, ns), tt = 0:(ns - 1)){
  m<-ncol(R)


  MM <- matrix(nrow = ncol(R), ncol = ns)
  mm <- matrix(nrow = ncol(R), ncol = ns)
  vv <- matrix(nrow = ncol(R), ncol = ns)
  time<-tt/max(tt)
  dt <- diff(time)

  Chol<-chol(R)

  temp<-MASS::mvrnorm(n =(ns), mu=rep(0,m), Sigma=((t(Chol)%*%Chol)*(dt[1])))

  for (i in 1:m){
    MM[i,c(1:ns)] <- cumsum(c(temp[,i]) )
    mm[i,c(1:ns)] <- MM[i,c(1:ns)] + rnorm(ns, 0, sqrt(vp/nn))
    vv[i,c(1:ns)]<-rep(vp,(ns))
  }

  MM<-MM+anc
  mm<-mm+anc
  
  List<-list()
  for (i in 1:m){
    List[[i]]<-paleoTS::as.paleoTS(mm = mm[i,], vv = vv[i,], nn = nn, tt = time, MM = MM[i,], label = "Created by sim.multi.BM", reset.time = FALSE)
  }

  if (m==2) yy<-make.multivar.evoTS(List[[1]], List[[2]])
  if (m==3) yy<-make.multivar.evoTS(List[[1]], List[[2]], List[[3]])
  if (m==4) yy<-make.multivar.evoTS(List[[1]], List[[2]], List[[3]], List[[4]])
  if (m==5) yy<-make.multivar.evoTS(List[[1]], List[[2]], List[[3]], List[[4]], List[[5]])


  return(yy)
}
