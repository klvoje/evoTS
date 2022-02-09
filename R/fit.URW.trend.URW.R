#' @title Fit three models (Unbiased Random Walk, General Random Walk (trend), Unbiased Random Walk) to three segments of an evolutionary sequence (time-series).
#'
#' @description Function to find maximum likelihood solutions to three models to three different segments of an evolutionary sequence.
#'
#' @param y an univariate evoTS object.
#'
#' @param minb the minimum number of samples within a segment to consider
#'
#' @param pool logical indicating whether to pool variances across samples
#'
#' @param silent if TRUE, less information is printed to the screen as the model is fit
#'
#' @param hess logical, indicating whether to calculate standard errors from the Hessian matrix.
#'
#' @param shift.point The samples that split the time-series into three segments. The samples are passed to the argument as a vector. Default is NULL, which means all possible switch points will be assessed constrained by how minib is defined.
#'#'
#'@return
#'\item{logL}{the log-likelihood of the optimal solution}
#'\item{AICc}{AIC with a correction for small sample sizes}
#'\item{parameters}{parameter estimates}
#'\item{modelName}{abbreviated model name}
#'\item{method}{Joint consideration of all samples}
#'\item{K}{number of parameters in the model}
#'\item{n}{the number of observations/samples}
#'\item{all.logl}{log-likelihoods for all tested partitions of the series into segments. Will return a single value if shift points have been given}
#'\item{GG}{matrix of indices of initial samples of each tested segment configuration; each column of GG corresponds to the elements of all.logl}
#'
#'@note The models have been implemented to be compatible with the joint parameterization routine in the package paleoTS. The optimization is therefore fit using the actual sample values, with the autocorrelation among samples accounted for in the log-likelihood function. The joint distribution of sample means is multivariate normal, with means and variance-covariances determined by evolutionary parameters and sampling errors.
#'
#'@author Kjetil Lysne Voje
#'
#'@references Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond. \emph{Paleobiology} 32:578–601.
#'
#'@export
#'
#'@examples
#'## Generate a evoTS objects by simulating a multivariate dataset
#'x <- paleoTS::sim.GRW(60)
#'
#'## Fit a the model to the data by defining shift points.
#'fit.URW.trend.URW(x, shift.point = c(20,40))
#'

fit.URW.trend.URW<-function (y, minb = 7, pool = TRUE, silent = FALSE, hess = FALSE, shift.point=NULL)
{
  ns <- length(y$mm)
  ng <- 3

  if(is.numeric(shift.point) == TRUE) GG <-shift.point else GG <- shifts(ns, ng, minb = minb)
  GG<-as.matrix(GG)
  if (ncol(GG) == 1) print("Fitting the model for a user-defined switchpoint") else print("Searching for all possible switchpoints in timeseries")

  #Define number of shift points:
  nc <- ncol(GG)
  if (!silent)
    cat("Total # hypotheses: ", nc, "\n")
  #Create empty list
  wl <- list()
  #create array with length = to switch points and where every entry ) -Inf
  logl <- array(-Inf, dim = nc)

  #start loop for estimating maximum likelihood parameters for each data set defined by the switch points
  for (i in 1:nc) {
    if (!silent)
      cat(i, " ")
    #defines which data point in the time series that belong to each of the two sets
    gg <- shift2gg(GG[, i], ns)
        w <- opt.joint.URW.trend.URW(y, gg, pool = pool, hess = hess)
    logl[i] <- w$logL
    wl[[i]] <- w
  }

  if (!silent)
    cat("\n")
  winner <- which.max(logl)
  ww <- wl[[winner]]
  ss <- GG[, winner]
  names(ss) <- paste("shift", 1:(ng - 1), sep = "")
  ww$parameters <- append(ww$parameters, ss)
  ww$all.logl <- logl
  ww$GG <- GG
  return(ww)
}
