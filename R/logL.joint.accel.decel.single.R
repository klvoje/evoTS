#' @title Log-likelihoods for evolutionary models
#'
#' @description Returns log-likelihood for a multivariate Unbiased Random Walk model with
#'   accelerating or decelerating rates of evolution through time.
#'
#' @param init.par initial (starting) parameters values
#'
#' @param y vector containing all trait values from all traits
#'
#' @param m number of traits
#'
#' @param n number of populations
#'
#' @param anc.values initial values for the ancestral trait values
#'
#' @param ta matrix of minimum times for each pair of time points
#'
#' @param se_vec pre-computed sampling error vector
#'
#' @details In general, users will not be access these functions directly, but instead use the
#'   optimization functions, which use these functions to find the best-supported parameter values.
#'
#'@return The log-likelihood of the parameter estimates, given the data.
#'
#'@author Kjetil Lysne Voje


### v.1.4 ###
# `yy` is replaced by `se_vec`, a sampling-error vector pre-computed once by the calling opt.* function.
logL.joint.accel.decel.single.R <- function(init.par, y, m, n, anc.values, ta, se_vec)
{
  m <- length(anc.values)
  chol <- diag(c(rep(0, m)))
  diag(chol) <- c(init.par[1:m])
  locations.R <- which(chol == 0, arr.ind = T)
  location.upper.tri.R <- which(locations.R[,1] < locations.R[,2])

  upper.first <- init.par[(m+1):(m+length(location.upper.tri.R))]

  for (i in 1:m){
    chol[locations.R[,1][location.upper.tri.R[i]], locations.R[,2][location.upper.tri.R[i]]] <- upper.first[i]
  }

  M.init <- init.par[(m+length(location.upper.tri.R)+1):(m+length(location.upper.tri.R)+m)]
### v.1.4 ###
# v.1.3 built M with a per-trait loop into M_temp, then transposed
# and flattened it. Replaced with a single vectorized call.
  M <- rep(M.init, each = n)  # vectorized: replaces the M_temp loop

### v.1.4 ###
# v.1.3 recomputed the full pairwise-min time matrix with outer()
# on every call (commented-out and live variants both present). Now built
# from the pre-computed `ta` argument instead.
  # C computed using pre-computed ta, avoiding recomputation of outer() on every call
  C <- (exp(init.par[length(init.par)] * ta) - 1) / init.par[length(init.par)]

  V  <- matrix(0, nrow = length(M), ncol = length(M))
  VV <- V + kronecker(t(chol) %*% chol, C)

### v.1.4 ###
# v.1.3 recomputed sample.var from yy$vv/yy$nn via a per-trait loop
# every call; now added directly from the pre-computed se_vec (see signature).
  diag(VV) <- diag(VV) + se_vec  # pre-computed sampling error

  y <- as.vector(y)
  S <- mvtnorm::dmvnorm(y, mean = M, sigma = VV, log = TRUE)
}
