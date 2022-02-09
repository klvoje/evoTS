#' @title Fit three Unbiased Random Walk models separated by punctuations (Bokma model) to an evolutionary sequence (time-series).
#'
#' @description Function to find maximum likelihood solutions to the Bokma model with two punctuations to an evolutionary sequence. Most users will interact with this function through the function fit.2.PE.Bokma.
#'
#' @param y an univariate evoTS object.
#'
#' @param gg numeric vector indicating membership of each sample in segments.
#'
#' @param pool logical indicating whether to pool variances across samples
#'
#' @param hess logical, indicating whether to calculate standard errors from the Hessian matrix.
#'
#' @param meth optimization method, passed to function optim. Default is "L-BFGS-B".
#'#'
#'@return
#'\item{logL}{the log-likelihood of the optimal solution}
#'\item{AICc}{AIC with a correction for small sample sizes}
#'\item{parameters}{parameter estimates}
#'\item{modelName}{abbreviated model name}
#'\item{method}{Joint consideration of all samples}
#'\item{K}{number of parameters in the model}
#'\item{n}{the number of observations/samples}
#'
#'@note The models have been implemented to be compatible with the joint parameterization routine in the package paleoTS. The optimization is therefore fit using the actual sample values, with the autocorrelation among samples accounted for in the log-likelihood function. The joint distribution of sample means is multivariate normal, with means and variance-covariances determined by evolutionary parameters and sampling errors.
#'
#'@author Kjetil Lysne Voje
#'
#'@references Bokma, F. 2002. Detection of punctuated equilibrium from molecular phylogenies. \emph{Journal of Evolutionary Biology} 15:1048–1056
#'
#'@export

opt.joint.2.Bokma<-function (y, gg, pool = TRUE, hess = FALSE, meth="L-BFGS-B")
{
  if (pool) {
    y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  }

  cl = list(fnscale = -1)

  small <- 1e-08
  y$tt <- y$tt - min(y$tt)
  p0 <- array(dim = 3)
  p0[1]<- y$mm[1]
  p0[2] <- paleoTS::mle.URW(paleoTS::sub.paleoTS(y, ok = gg == 1))
  if (p0[2] <= small) p0[2] <- 100 * small
  p0[3] <- (paleoTS::sub.paleoTS(y, ok = gg == 2))$mm[1]
  p0[4] <- paleoTS::mle.URW(paleoTS::sub.paleoTS(y, ok = gg == 2))
  if (p0[4] <= small) p0[4] <- 100 * small
  p0[5] <- (paleoTS::sub.paleoTS(y, ok = gg == 3))$mm[1]
  p0[6] <- paleoTS::mle.URW(paleoTS::sub.paleoTS(y, ok = gg == 3))
  names(p0) <- c("anc", "vstep_1", "opt.2", "vstep_2", "opt.3", "vstep_3")
  K <- 8

  cl$ndeps <- p0/100

      w <- optim(p0, fn = logL.joint.2.Bokma, gg = gg, method = meth,
                 lower = c(NA, 0, NA, 0), control = cl,
                 hessian = hess, y = y)

  if (hess)
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- paleoTS::as.paleoTSfit(logL = w$value, parameters = w$par, modelName = paste("Bokma model with two punctuations",
                                                                            sep = "-"), method = "Joint", K = K, n = length(y$mm),
                      se = w$se)
  return(wc)
}
