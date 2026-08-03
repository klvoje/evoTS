#' @title Log-likelihoods for evolutionary models
#'
#' @description Returns log-likelihood for a multivariate Ornstein-Uhlenbeck model
#'   with user-defined A and R matrices.
#'
#' @param init.par initial (starting) parameters values
#'
#' @param yy a multivariate evoTS object
#'
#' @param A.user the pull matrix.
#'
#' @param R.user the drift matrix.
#'
#' @param locations.A location (row and column) of parameters (elements) in the A matrix that is estimated
#'
#' @param location.diag.A location (row and column) of parameters (elements) in the diagonal of the A matrix that is estimated
#'
#' @param location.upper.tri.A location (row and column) of parameters (elements) in the upper triangle of the A matrix that is estimated
#'
#' @param location.lower.tri.A location (row and column) of parameters (elements) in the lower triangle of the A matrix that is estimated
#'
#' @param locations.R location (row and column) of parameters (elements) in the R matrix that is estimated
#'
#' @param location.diag.R location (row and column) of parameters (elements) in the diagonal of the R matrix that is estimated
#'
#' @param location.upper.tri.R location (row and column) of parameters (elements) in the upper triangle of the R matrix that is estimated
#'
#' @param ta matrix of minimum times for each pair of time points
#'
#' @param tij matrix of absolute time differences for each pair of time points
#'
#' @param time_vec rescaled time vector
#'
#' @param se_vec pre-computed sampling error vector
#'
#' @details In general, users will not access this function directly, but instead use the
#'   optimization functions, which use this function to find the best-supported parameter
#'   values.
#'
#'@return The log-likelihood of the parameter estimates, given the data.
#'
#'@author Kjetil Lysne Voje


### new in v.1.4 ###
# Signature changed: `ta` (pairwise-minimum time matrix), `tij` (pairwise
# absolute time-difference matrix), `time_vec` (rescaled time vector), and
# `se_vec` (pre-computed sampling error vector) are new arguments, computed
# once by the calling fit.multivariate.OU.user.defined function instead of
# being rebuilt from `yy` on every call. `yy` itself is retained, but is now
# only used for m/X/y below.
logL.joint.multi.OUOU.user <- function(init.par, yy, A.user, R.user,
                                           locations.A, location.diag.A,
                                           location.upper.tri.A, location.lower.tri.A,
                                           locations.R, location.diag.R, location.upper.tri.R,
                                           ta, tij, time_vec, se_vec) {

  m   <- ncol(yy$xx)
  X   <- yy$xx
  y   <- as.matrix(as.vector(X))
  n_t <- length(time_vec)

  ### -----------------------------------------------------------------------
  ### Parameter unpacking into A and Chol matrices
  ### (unchanged logic from original logL.joint.multi.OUOU.user)
  ### -----------------------------------------------------------------------

  A <- diag(rep(0.00000000001, m))

  for (i in 1:length(location.diag.A)) {
    A[locations.A[location.diag.A][i], locations.A[location.diag.A][i]] <- init.par[i]
  }

  if (pracma::isempty(location.upper.tri.A) == FALSE) {
    for (i in 1:length(location.upper.tri.A)) {
      A[locations.A[, 1][location.upper.tri.A][i], locations.A[, 2][location.upper.tri.A][i]] <- init.par[(length(location.diag.A) + i)]
    }
  } else location.upper.tri.A <- NULL

  if (pracma::isempty(location.lower.tri.A) == FALSE) {
    for (i in 1:length(location.lower.tri.A)) {
      A[locations.A[, 1][location.lower.tri.A][i], locations.A[, 2][location.lower.tri.A][i]] <- init.par[(length(location.diag.A) + length(location.upper.tri.A) + i)]
    }
  } else location.lower.tri.A <- NULL

### v.1.4 ###
# New: eigendecomposition computed once and reused (v.1.3 called
# eigen(A) a second time later to get D, repeating the decomposition).
  # Eigendecomposition of A: computed once here (v.1.3 called eigen(A) twice)
  eig_A <- eigen(A)
### end v.1.4 ###
  P     <- eig_A$vectors
  D     <- diag(eig_A$values)

  Chol <- diag(rep(0, m))
  for (i in 1:length(location.diag.R)) {
    Chol[locations.R[location.diag.R][i], locations.R[location.diag.R][i]] <- init.par[(length(location.diag.A) + length(location.upper.tri.A) + length(location.lower.tri.A) + i)]
  }

  if (pracma::isempty(location.upper.tri.R) == FALSE) {
    for (i in 1:length(location.upper.tri.R)) {
      Chol[locations.R[, 1][location.upper.tri.R][i], locations.R[, 2][location.upper.tri.R][i]] <- init.par[(length(location.diag.A) + length(location.upper.tri.A) + length(location.lower.tri.A) + length(location.diag.R) + i)]
    }
  } else location.upper.tri.R <- NULL

  ### Theta (optimal trait values) ###
  optima <- c(init.par[(length(location.diag.A) + length(location.upper.tri.A) + length(location.lower.tri.A) + length(location.diag.R) + length(location.upper.tri.R) + 1):(length(location.diag.A) + length(location.upper.tri.A) + length(location.lower.tri.A) + length(location.diag.R) + length(location.upper.tri.R) + m)])

  ### The ancestral trait values ###
  anc <- c(init.par[(length(location.diag.A) + length(location.upper.tri.A) + length(location.lower.tri.A) + length(location.diag.R) + length(location.upper.tri.R) + m + 1):(length(location.diag.A) + length(location.upper.tri.A) + length(location.lower.tri.A) + length(location.diag.R) + length(location.upper.tri.R) + m + m)])

### v.1.4 ###
# To increase speed: v.1.3  computed the time vector, the per-time-point exp(-A*t)
# array, M, and VV3 by rebuilding everything from `yy$tt` inside this function
# on every call, using explicit i/j/k/l loops throughout (including a nested
# k,l loop for the VCV integral and a tmp.VV intermediate array for block
# assembly). This has been rewritten to use the pre-computed ta/tij/time_vec/se_vec
# arguments, vectorised matrix expressions in place of the innermost loops,
# and to add tryCatch/finite-value guards (see end of function) that were not
# present in v.1.3.

### -----------------------------------------------------------------------
### Optimization: pre-compute constants that are reused across time pairs
### -----------------------------------------------------------------------

  # P_inv: computed once instead of inside the i,j loop
  P_inv <- solve(P)

  # d: eigenvalues of A as a plain vector
  d <- diag(D)

  # sum_d: m x m matrix of all pairwise eigenvalue sums (dk + dl).
  # Replaces the innermost k,l double loop in the original.
  sum_d <- outer(d, d, "+")

  # right.side: the P^{-1} R P^{-T} factor in the integral; constant across all
  # time pairs so it is computed once.
  right.side <- P_inv %*% (Chol %*% t(Chol)) %*% t(P_inv)

  ### -----------------------------------------------------------------------
  ### Calculate expected trait values M
  ### -----------------------------------------------------------------------

  anc_minus_theta <- anc - optima
  M <- matrix(NA, nrow = m, ncol = n_t)

  for (i in 1:n_t) {
    eAt_i  <- P %*% diag(exp(-d * time_vec[i])) %*% P_inv
    M[, i] <- eAt_i %*% anc_minus_theta + optima
  }
  M <- c(t(M))   # vectorize row-by-row to match ordering of y

  ### -----------------------------------------------------------------------
  ### Compute variance-covariance matrix VV3(eq. 8 and 9 from Suppl. of Clavel et al. 2015)
  ###
  ### v.1.4 vs. v.1.3:
  ###   - No tmp.VV intermediate array; VV3 is filled directly.
  ###   - left.side computed by outer() instead of a k,l double loop.
  ###   - exp decay vector computed by a single vectorized exp() call
  ###     (also fixes variable-shadowing bug: original used 'm' as loop
  ###     variable, overwriting the trait-count variable).
  ###   - row/col block offsets pre-computed outside the inner loop.
  ### -----------------------------------------------------------------------

  VV3      <- matrix(0, nrow = n_t * m, ncol = n_t * m)
  row_base <- (0:(m - 1)) * n_t
  col_base <- row_base

  for (a in 1:n_t) {
    for (b in 1:n_t) {

      # Vectorized left.side: replaces the original double loop over k and l.
      left.side  <- (1 - exp(-sum_d * ta[a, b])) / sum_d
      left.right <- left.side * right.side
      integ      <- P %*% left.right %*% t(P)

      # Vectorized exp decay: replaces the original inner loop over traits.
      # Variable name 'exp_diag_b' avoids shadowing the outer 'm'.
      exp_diag_b <- exp(-d * tij[a, b])
      exp_mat_T  <- t(P %*% diag(exp_diag_b) %*% P_inv)

      tmp <- integ %*% exp_mat_T

      # Direct placement into VV3 (replaces tmp.VV array and assembly step)
      row_idx <- row_base + a
      col_idx <- col_base + b
      VV3[row_idx, col_idx] <- t(tmp)
    }
  }

  ### Add pre-computed sampling error to the diagonal ###
  diag(VV3) <- diag(VV3) + se_vec
  VV3 <- (VV3 + t(VV3)) / 2

  if (any(!is.finite(M)) || any(!is.finite(VV3))) return(-1e20)

  S <- tryCatch(
    mvtnorm::dmvnorm(t(y), mean = M, sigma = VV3, log = TRUE),
    error = function(e) -1e20
  )
  if (!is.finite(S)) return(-1e20)
  S
}
