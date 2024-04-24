#' @title Fit an Ornstein-Uhlenbeck model with an optimum that evolves according to a Unbiased Random Walk.
#'
#' @description Function to find maximum likelihood solutions to an Ornstein-Uhlenbeck model with an optimum that evolves according to a Unbiased Random Walk.
#'
#' @param y an univariate paleoTS object.
#'
#' @param pool logical indicating whether to pool variances across samples
#'
#' @param meth optimization method, passed to function optim. Default is "L-BFGS-B".
#'
#' @param hess logical, indicating whether to calculate standard errors from the Hessian matrix.
#'
#' @param iterations the number of times the optimization method is run from different starting points. Default is NULL, meaning the optimization is run once.
#'
#' @param iter.sd defines the standard deviation of the Gaussian distribution from which starting values for the optimization routine is run. Default is 1.
#'
#' @param opt.anc logical, indicating whether the the ancestral trait state is at the optimum.
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
#'@references Hansen, T. F., Pienaar, J. & Orzack, S. H. 2008. A Comparative Method for Studying Adaptation to a Randomly Evolving Environment. \emph{Evolution} 62:1965–1977.
#'
#'@export
#'
#'@examples
#'## Generate a paleoTS object by simulating a univariate evolutionary sequence
#'x <- paleoTS::sim.GRW(60)
#'
#'## Fit the model
#'opt.joint.OUBM(x)
#'

opt.joint.OUBM<-function (y, pool = TRUE, meth = "L-BFGS-B", hess = FALSE, iterations=NULL, iter.sd=NULL, opt.anc=TRUE)
{
  if (pool)
    y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  cl <- list(fnscale = -1)
  if (y$tt[1] != 0)
    stop("Initial time must be 0.  Use as.paleoTS() or read.paleoTS() to correctly process ages.")
  
  if (is.numeric(iterations)) cat("The optimization method is executed from multiple different starting points. Number of iterations:", iterations)

  w0 <- paleoTS::mle.GRW(y)
  bm0<-paleoTS::mle.GRW(y)
  halft <- (y$tt[length(y$tt)] - y$tt[1])/4
  if (opt.anc == TRUE) p0 <- c(y$mm[1], w0[2]/10, log(2)/halft, bm0[2]) else p0 <- c(y$mm[1], w0[2]/10, y$mm[1], log(2)/halft, bm0[2])

  if (opt.anc == TRUE) names(p0) <- c("anc/theta.0", "vstep.trait", "alpha", "vstep.opt") else names(p0) <- c("anc", "vstep.trait", "theta.0", "alpha", "vstep.opt")

  base.init.par<-p0
  if (is.numeric(iterations)) {
    if(is.numeric(iter.sd) == FALSE) iter.sd <-1
    ww<-rep(NA, iterations)

    init.par<-matrix(NA, nrow=iterations, ncol=length(base.init.par))
    for (i in 1: iterations){
        init.par[i,]<-c(rnorm(length(base.init.par), base.init.par, iter.sd))

        if (opt.anc == TRUE) init.par[i,c(2,3,4)]<-abs(init.par[i,c(2,3,4)]) else init.par[i,c(2,4,5)]<-abs(init.par[i,c(2,4,5)])
        if (opt.anc == TRUE) {
          tmp <- optim(init.par[i,], fn = logL.joint.OU.BM, control = cl, method = "L-BFGS-B",
                     lower = c(NA, 1e-08, 1e-08, 1e-08), hessian = hess,
                     y = y, opt.anc = opt.anc)}
        if (opt.anc == FALSE) {
          tmp <- optim(init.par[i,], fn = logL.joint.OU.BM, control = cl, method = "L-BFGS-B",
                       lower = c(NA, 1e-08, NA, 1e-08, 1e-08), hessian = hess,
                       y = y, opt.anc = opt.anc)}
        ww[i]<-tmp$value
    }

    best.run<-match(max(na.exclude(ww)), ww)
    init.par<-init.par[best.run,]
    iter<-length(which(!is.na(ww)))
  }

  if (is.numeric(iterations)) p0<-init.par else p0<-p0
  if (opt.anc == TRUE) names(p0) <- c("anc/theta.0", "vstep.trait", "alpha", "vstep.opt") else names(p0) <- c("anc", "vstep.trait", "theta.0", "alpha", "vstep.opt")

  if (is.null(cl$ndeps))
    cl$ndeps <- abs(p0/10000)
  cl$ndeps[cl$ndeps == 0] <- 1e-08
  if (meth == "L-BFGS-B")
    if (opt.anc == TRUE)
      w <- optim(p0, fn = logL.joint.OU.BM, control = cl, method = "L-BFGS-B",
                 lower = c(NA, 1e-08, 1e-08, 1e-08), hessian = hess,
                 y = y, opt.anc = opt.anc)
  else
   w <- optim(p0, fn = logL.joint.OU.BM, control = cl, method = "L-BFGS-B",
               lower = c(NA, 1e-08, NA, 1e-08, 1e-08), hessian = hess,
               y = y, opt.anc = opt.anc)
  else w <- optim(p0, fn = logL.joint.OU.BM, control = cl, method = meth,
                  hessian = hess, y = y, opt.anc = opt.anc)

  if (is.numeric(iterations) == FALSE) {
    iter<-NA
  }

  if (opt.anc == TRUE) K<-4 else K<-5
  if (opt.anc == TRUE) modelname<-"OU model with moving optimum (ancestral state at optimum)" else modelname<-"OU model with moving optimum"


  if (hess)
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- as.evoTSfit.OUBM(logL = w$value, parameters = w$par, modelName = modelname,
                      method = "Joint", K = K, n = length(y$mm), se = w$se)
  return(wc)
}
