#' @title Simulate an Ornstein-Uhlenbeck process with optimum changing according to an unbiased random walk
#'
#' @description Function to simulate an Ornstein-Uhlenbeck evolutionary sequence data set with an optimum moving according to an unbiased random walk.
#'
#' @param ns number of samples in time-series
#'
#' @param anc the ancestral trait values
#'
#' @param theta.0 the ancestral value for the optimum
#'
#' @param alpha strength of attraction to the optimum
#'
#' @param vstep.trait step variance of the trait
#'
#' @param vstep.opt step variance of the optimum
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
#'##Simulate data
#'x<-sim.OUBM(50, theta.0 = 5, alpha = 0.6, vstep.opt = 0.5)
#'
#'## plot the data
#'plot(x)


sim.OUBM<-function (ns = 20, anc = 0, theta.0 = 0, alpha = 0.3, vstep.trait = 0.1, vstep.opt = 0.1,
          vp = 1, nn = rep(20, ns), tt = 0:(ns - 1))
{
  vs<-vstep.trait
  vo<-vstep.opt

  ou.MM <- function(anc, theta.0, aa, tt) theta.0*(1 - exp(-aa*tt)) + anc*exp(-aa*tt)
  ou.VV <-function(vo, vs, aa, tt)  ((vo+vs)/(2 * aa)) * (1 - exp(-2 * aa * tt)) + (vo*tt*(1-(2*(1-exp(-aa*tt)))))/(aa*tt)

  MM <- array(dim = ns)
  mm <- array(dim = ns)
  vv <- array(dim = ns)
  dt <- diff(tt)
  MM[1] <- anc
  x <- stats::rnorm(nn[1], mean = MM[1], sd = sqrt(vp))
  mm[1] <- mean(x)
  vv[1] <- var(x)
  for (i in 2:ns) {
    ex <- ou.MM(MM[i - 1], theta.0, alpha, dt[i - 1])
    vx <- ou.VV(vo, vs, alpha, dt[i - 1])
    MM[i] <- rnorm(1, ex, sqrt(vx))
    x <- stats::rnorm(nn[i], mean = MM[i], sd = sqrt(vp))
    mm[i] <- mean(x)
    vv[i] <- var(x)
  }
  gp <- c(anc, theta.0, alpha, vstep.trait, vstep.opt)
  names(gp) <- c("anc", "theta.0", "alpha", "vstep.trait", "vstep.opt")
  res <- paleoTS::as.paleoTS(mm = as.vector(mm), vv = as.vector(vv),
                    nn = nn, tt = tt, MM = MM, genpars = gp, label = "Created by sim.OUBM()",
                    reset.time = FALSE)
  return(res)
}
