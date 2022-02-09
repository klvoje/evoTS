#' @title Fit two Unbiased Random Walk models separated by punctuations (Bokma model) to an evolutionary sequence (time-series).
#'
#' @description Function to find maximum likelihood solutions to the Bokma model with a single punctuation to an evolutionary sequence.
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
#'@references Bokma, F. 2002 Detection of punctuated equilibrium from molecular phylogenies. \emph{J Evolution Biol} 15:1048–1056.
#'
#'@export
#'
#'@examples
#'## Generate a evoTS objects by simulating a multivariate dataset
#'x <- paleoTS::sim.GRW(60)
#'
#'## Fit a the model to the data by defining shift points.
#'fit.Bokma(x, shift.point = 20)
#'

fit.Bokma<-function (y, minb = 10, pool = TRUE, silent = FALSE, hess = FALSE, shift.point=NULL)
{
  ng <- 2
  ns <- length(y$mm)

  if(is.numeric(shift.point) == TRUE) GG <-shift.point else GG <- shifts(ns, ng, minb = minb)
  GG<-as.matrix(GG)
  if (ncol(GG) == 1) print("Fitting the model for a user-defined switchpoint") else print("Searching for all possible switchpoints in timeseries")

  nc <- ncol(GG)
  if (!silent)
    cat("Total # hypotheses: ", nc, "\n")
  wl <- list()
  logl <- array(-Inf, dim = nc)
  for (i in 1:nc) {
    if (!silent)
      cat(i, " ")
    gg <- shift2gg(GG[, i], ns)
    w <- opt.joint.Bokma(y, gg , pool = pool, hess = hess)
    logl[i] <- w$logL
    wl[[i]] <- w
  }
  cat("\n")
  winner <- which.max(logl)
  ww <- wl[[winner]]
  ss <- GG[, winner]
  names(ss) <- paste("shift", 1, sep = "")
  ww$parameters <- append(ww$parameters, ss)
  ww$all.logl <- logl
  ww$GG <- GG
  return(ww)
}
