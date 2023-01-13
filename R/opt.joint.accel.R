#' @title Fit an Unbiased Random Walk with an accelerating rate of change through time.
#'
#' @description Function to find maximum likelihood solutions to a Unbiased Random Walk with an accelerating or decelerating rate of change through time.
#'
#' @param y an univariate evoTS object.
#'
#' @param pool logical indicating whether to pool variances across samples
#'
#' @param meth optimization method, passed to function optim. Default is "L-BFGS-B".
#'
#' @param hess logical, indicating whether to calculate standard errors from the Hessian matrix.
#'
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
#'@export
#'
#'@examples
#'## Generate a paleoTS object by simulating a univariate evolutionary sequence
#'y <- paleoTS::sim.GRW(30)
#'
#'## Fit the model
#'opt.joint.accel(y)
#'
opt.joint.accel<-function (y, pool = TRUE, meth = "L-BFGS-B", hess = FALSE)
{
  if (pool)
    y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)

  cl = list(fnscale = -1)

  if (y$tt[1] != 0)
    stop("Initial time must be 0.  Use as.paleoTS() or read.paleoTS() to correctly process ages.")
  p0 <- array(dim = 3)
  p0[1] <- y$mm[1]
  p0[2] <- min(c(paleoTS::mle.URW(y), 1e-07))
  p0[3] <- 1
  names(p0) <- c("anc", "vstep", "r")
  if (is.null(cl$ndeps))
    cl$ndeps <- abs(p0/10000)
  cl$ndeps[cl$ndeps == 0] <- 1e-08
  if (meth == "L-BFGS-B")
    w <- optim(p0, fn = logL.joint.accel_decel, control = cl, method = meth,
               lower = c(NA, 0, +0.000001), upper = (c(NA, NA, NA)), hessian = hess, y = y)
  else w <- optim(p0, fn = logL.joint.accel_decel, control = cl, method = meth, lower = c(NA, 0, NA),
                  upper = (c(NA, NA,0)), hessian = hess, y = y)

  if (hess)
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- paleoTS::as.paleoTSfit(logL = w$value, parameters = w$par, modelName = "Accel",
                      method = "Joint", K = 3, n = length(y$mm), se = w$se)
  return(wc)
}
