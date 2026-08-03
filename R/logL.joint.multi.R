#' @title Log-likelihoods for evolutionary models
#'
#' @description Returns log-likelihood for a multivariate Unbiased Random Walk model fitted
#'   to separate segments of a multivariate time series.
#'
#' @param init.par initial (starting) parameters values
#'
#' @param C distance matrix (list with two elements: C.1 and C.2)
#'
#' @param y vector containing all trait values from all traits
#'
#' @param m number of traits
#'
#' @param n number of populations
#'
#' @param anc.values initial values for the ancestral trait values
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
logL.joint.multi.R <- function(init.par, C, y, m, n, anc.values, se_vec)
{
  m <- length(anc.values)
  var.covar.1 <- diag(c(rep(0, m)))
  var.covar.2 <- diag(c(rep(0, m)))
  diag(var.covar.1) <- c(init.par[1:m])
  diag(var.covar.2) <- c(init.par[(m+1):(m+m)])
  locations.R <- which(var.covar.1 == 0, arr.ind = T)
  location.upper.tri.R <- which(locations.R[,1] < locations.R[,2])

  upper.first  <- init.par[(m+m+1):(m+m+length(location.upper.tri.R))]
  upper.second <- init.par[(m+m+length(location.upper.tri.R)+1):(m+m+length(location.upper.tri.R)+length(location.upper.tri.R))]

  for (i in 1:m){
    var.covar.1[locations.R[,1][location.upper.tri.R[i]], locations.R[,2][location.upper.tri.R[i]]] <- upper.first[i]
    var.covar.2[locations.R[,1][location.upper.tri.R[i]], locations.R[,2][location.upper.tri.R[i]]] <- upper.second[i]
  }

  cholesky <- array(c(var.covar.1, var.covar.2), dim = c(m, m, 2))

  M.init <- init.par[(m+m+length(location.upper.tri.R)+length(location.upper.tri.R)+1):(m+m+length(location.upper.tri.R)+length(location.upper.tri.R)+m)]
### v.1.4 ###
# V.1.3 built M with a per-trait loop into M_temp, then transposed
# and flattened it. Replaced with a single vectorized call.
  M <- rep(M.init, each = n)  # vectorized: replaces the M_temp loop

  V     <- matrix(0, nrow = length(M), ncol = length(M))
  tmp.1 <- V + kronecker(t(cholesky[,,1]) %*% cholesky[,,1], C[[1]])
  tmp.2 <- V + kronecker(t(cholesky[,,2]) %*% cholesky[,,2], C[[2]])
  VV    <- tmp.1 + tmp.2

### v.1.4 ###
# V.1.3 recomputed sample.var from yy$vv/yy$nn via a per-trait loop
# every call; now added directly from the pre-computed se_vec.
  diag(VV) <- diag(VV) + se_vec  # pre-computed sampling error

### v.1.4 ###
# Guards against non-finite mean/covariance (e.g. during bad optimizer
# proposals) by returning a large penalty instead of letting dmvnorm error out.
  if (any(!is.finite(M)) || any(!is.finite(VV))) return(-1e20)

  y <- as.vector(y)
### v.1.4 ###
# V.1.3 called dmvnorm directly and returned its result unguarded.
# Now wrapped in tryCatch and finite-checked, both falling back to a large
# penalty (-1e20) rather than crashing/propagating NA into the optimizer.
  S <- tryCatch(
    mvtnorm::dmvnorm(y, mean = M, sigma = VV, log = TRUE),
    error = function(e) -1e20
  )
  if (!is.finite(S)) return(-1e20)
  S
}
