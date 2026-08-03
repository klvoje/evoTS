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

  MM <- matrix(nrow = m, ncol = ns)
  mm <- matrix(nrow = m, ncol = ns)
  vv <- matrix(nrow = m, ncol = ns)
  time <- tt / max(tt)
  dt   <- diff(time)

### v.1.4 ###
# Bug fix. Version 1.3 drew all ns increments from a
# single mvrnorm() call using Sigma = R * dt[1] (only the first inter-sample
# interval), so every increment was scaled identically regardless of the
# actual (possibly unequal) spacing of later time points,

# Also, now each of the ns-1 increments is drawn with its own dt[s]-scaled
# covariance (correct for unevenly-spaced time series), and the first time
# point is fixed exactly at the ancestral value.
  
  # Start at ancestral values; simulate ns-1 increments, each scaled by its own dt
  MM[, 1] <- 0
  for (s in 1:(ns - 1)) {
    MM[, s + 1] <- MASS::mvrnorm(1, mu = rep(0, m), Sigma = R * dt[s])
  }
  for (i in 1:m) {
    MM[i, ] <- cumsum(MM[i, ])
  }
  MM <- MM + anc

  for (i in 1:m) {
    mm[i, ] <- MM[i, ] + rnorm(ns, 0, sqrt(vp / nn))
    vv[i, ] <- rep(vp, ns)
  }

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
