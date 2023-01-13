#' @title Simulate an Unbiased Random Walk with an accelerating or decelerating rate of change through time.
#'
#' @description Function to simulate an evolutionary sequence data set according to an Unbiased Random Walk with an accelerating or decelerating rate of change through time.
#'
#' @param ns number of samples in time-series
#'
#' @param vs step variance of the trait
#'
#' @param r the parameter controlling the exponential decay (if negative) or increase (if positive) of the rate (vs) through time.
#'
#' @param vp phenotypic variance of each sample
#'
#' @param nn 	vector of the number of individuals in each sample (identical sample sizes for all time-series is assumed)
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
#'x<-sim.accel.decel(40, r=-0.5)
#'
#'## Plot the data
#'plotevoTS(x)

sim.accel.decel<-function (ns = 20, vs = 0.5, r = 0.2, vp = 0.2, nn = rep(20, ns),
          tt = 0:(ns - 1))
{
  MM <- array(dim = ns)
  mm <- array(dim = ns)
  vv <- array(dim = ns)
  dt <- diff(tt)
  inc <- stats::rnorm(ns - 1, 0, sqrt(vs * exp(r*tt)))
  MM <- cumsum(c(0, inc))
  mm <- MM + rnorm(ns, 0, sqrt(vp/nn))
  vv <- rep(vp, ns)
  gp <- c(vs, r)
  names(gp) <- c("vstep", "r")
  res <- paleoTS::as.paleoTS(mm = mm, vv = vv, nn = nn, tt = tt, MM = MM,
                    genpars = gp, label = "Created by sim.accel", reset.time = FALSE)
  return(res)
}
