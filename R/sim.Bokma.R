#' @title Simulate the Bokma model
#'
#' @description Function to simulate an evolutionary sequence data set with two Unbiased Random Walks separated by a punctuation event (jump in phenotype space)
#'
#' @param ns number of samples in time-series
#'
#' @param vs step variance of the trait
#'
#' @param vp phenotypic variance of each sample
#'
#' @param nn 	vector of the number of individuals in each sample (identical sample sizes for all time-series is assumed)
#'
#' @param anc vector of length two with the starting values for each of the two unbiased random walks.
#'
#' @param tt 	vector of sample times (ages
#'
#'@return An evolutionary sequence (time-series) data set (a paleoTS object)
#'
#'@author Kjetil Lysne Voje
#'
#'@export
#'
#'@examples
#'##Simulate an unbiased random walk where the rate decelerates through time.
#'x<-sim.Bokma(c(20,20))
#'
#'## Plot the data
#'plot(x)
sim.Bokma<-function (ns = c(20, 20), vs = c(0.1, 0.1), vp = 1,
                     nn = rep(30, sum(ns)), anc = c(0, 2), tt = 0:(sum(ns) - 1))

{
  MM.1 <- array(dim = ns[1])
  mm.1 <- array(dim = ns[1])
  vv.1 <- array(dim = ns[1])
  dt <- diff(tt[c(1:ns[1])])
  inc <- stats::rnorm(ns[1]-1, 0, sqrt(vs[1] * dt))
  MM.1 <- cumsum(c(0, inc))
  MM.1<-MM.1+anc[1]
  mm.1 <- MM.1 + stats::rnorm(ns[1], 0, sqrt(vp/nn[c(1:ns[1])]))
  vv.1 <- rep(vp, ns[1])

  MM.2 <- array(dim = ns[2])
  mm.2 <- array(dim = ns[2])
  vv.2 <- array(dim = ns[2])
  dt <- diff(tt[(ns[1]+1):(ns[2]*2)])
  inc <- stats::rnorm(ns[2]-1, 0, sqrt(vs[2] * dt))
  MM.2 <- cumsum(c(0, inc))
  MM.2<-MM.2+anc[2]
  mm.2 <- MM.2 + stats::rnorm(ns[2], 0, sqrt(vp/nn[(ns[1]+1):(ns[2]*2)]))
  vv.2 <- rep(vp, ns[2])


  gp <- c(vs[1], anc[1]-anc[2], vs[1])
  names(gp) <- c("vstep.1", "punctuation", "vstep.2")
  res <- paleoTS::as.paleoTS(mm = c(mm.1, mm.2), vv = c(vv.1, vv.2), nn = nn, tt = tt, MM = c(MM.1, MM.2),
                    genpars = gp, label = "Created by sim.Bokma", reset.time = FALSE)
  return(res)
}

