#' @title Fit three OU models to three segments of an evolutionary sequence (time-series).
#'
#' @description Function to find maximum likelihood solutions to three OU models to three different segments of an evolutionary sequence.
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
#' @param shift.point The samples that split the time-series into three segments. The samples are passed to the argument as a vector. Default is NULL, which means all possible shift points will be assessed constrained by how minb is defined.
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
#'@references Hunt, G., Bell, M. A. & Travis, M. P. Evolution toward a new adaptive optimum: Phenotypic evolution in a fossil stickleback lineage. \emph{Evolution} 62:700–710.
#'
#'@export
#'
#'@examples
#'## Generate a evoTS objects by simulating a multivariate dataset
#'x <- paleoTS::sim.GRW(60)
#'
#'## Fit a the model to the data by defining shift points.
#'fit.OU.OU.OU(x, shift.point = c(20,40))
#'

fit.OU.OU.OU<-function (y, minb = 10, pool = TRUE, silent = FALSE, hess = FALSE, shift.point=NULL)
{
  ng <- 3
  ns <- length(y$mm)

  if(is.numeric(shift.point) == TRUE) GG <-shift.point else GG <- shifts(ns, ng, minb = minb)
  GG<-as.matrix(GG)
  if (ncol(GG) == 1) print("Fitting the model for a user-defined shift point") else print("Searching for all possible shift points in timeseries")

  nc <- ncol(GG)
  if (!silent)
    cat("Total # hypotheses: ", nc, "\n")
  wl <- list()
  logl <- array(-Inf, dim = nc)
  for (i in 1:nc) {
    if (!silent)
      cat(i, " ")
    gg <- shift2gg(GG[, i], ns)
    w <- opt.joint.OU.OU.OU(y, gg , pool = pool, hess = hess)
    logl[i] <- w$logL
    wl[[i]] <- w
  }
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
