#' @title Log-likelihoods for evolutionary models
#'
#' @description Returns log-likelihood for a multivariate Ornstein-Uhlenbeck model.
#'
#' @param init.par initial (starting) parameters values
#'
#' @param yy a multivariate evoTS object
#'
#' @param A.matrix the pull matrix.
#'
#' @param R.matrix the drift matrix.
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
# `ta` (pairwise-minimum time matrix), `tij` (pairwise absolute time-difference matrix), `time_vec` (rescaled time vector), and
# `se_vec` (pre-computed sampling error vector) are new arguments, computed once by the calling fit.multivariate.OU* function instead of being rebuilt
# from `yy` on every call. `yy` itself is retained, but is now only used for m/X/y below.
logL.joint.multi.OUOU <- function(init.par, yy, A.matrix, R.matrix,
                                      ta, tij, time_vec, se_vec) {

  m   <- ncol(yy$xx)          # number of traits
  X   <- yy$xx
  y   <- as.matrix(as.vector(X))
  n_t <- length(time_vec)     # number of time points

  ### -----------------------------------------------------------------------
  ### Parameter unpacking (unchanged from original logL.joint.multi.OUOU)
  ### -----------------------------------------------------------------------

  if (A.matrix == "diag" & R.matrix == "diag") {

    A    <- diag(c(init.par[1:m]))
    P    <- eigen(A)$vectors
    D    <- diag(eigen(A)$values)
    Chol <- diag(init.par[(m + 1):(m * 2)])

    optima <- c(init.par[((m * 2) + 1):(m * 3)])
    anc    <- c(init.par[((m * 3) + 1):(m * 4)])
  }

  if (A.matrix == "diag" & R.matrix == "symmetric") {

    A    <- diag(c(init.par[1:m]))
    P    <- eigen(A)$vectors
    D    <- diag(eigen(A)$values)
    Chol <- diag(init.par[(m + 1):(m * 2)])

    nr.off.diag <- upper.tri(Chol)
    l.upp.tri   <- length(nr.off.diag[nr.off.diag == TRUE])
    Chol[upper.tri(Chol)] <- init.par[((length(diag(A)) * 2) + 1):((length(diag(A)) * 2) + l.upp.tri)]

    optima <- c(init.par[((length(diag(A)) * 2) + l.upp.tri + 1):((length(diag(A)) * 2) + l.upp.tri + m)])
    anc    <- c(init.par[((length(diag(A)) * 2) + l.upp.tri + 1 + m):((length(diag(A)) * 2) + l.upp.tri + m + m)])
  }

  if (A.matrix == "full" & R.matrix == "diag") {

    A <- diag(c(init.par[1:m]))
    nr.off.diag <- upper.tri(A)
    l.upp.tri   <- length(nr.off.diag[nr.off.diag == TRUE])
    A[upper.tri(A)] <- init.par[(length(diag(A)) + 1):(length(diag(A)) + l.upp.tri)]
    A[lower.tri(A)] <- init.par[(length(diag(A)) + l.upp.tri + 1):(length(diag(A)) + l.upp.tri + l.upp.tri)]
    P <- eigen(A)$vectors
    D <- diag(eigen(A)$values)

    Chol <- diag(init.par[(length(diag(A)) + l.upp.tri + l.upp.tri + 1):(length(diag(A)) + l.upp.tri + l.upp.tri + length(diag(A)))])

    optima <- c(init.par[(length(diag(A)) + l.upp.tri + l.upp.tri + length(diag(A)) + 1):(length(diag(A)) + l.upp.tri + l.upp.tri + length(diag(A)) + length(diag(A)))])
    anc    <- c(init.par[(length(diag(A)) + l.upp.tri + l.upp.tri + length(diag(A)) + length(diag(A)) + 1):(length(diag(A)) + l.upp.tri + l.upp.tri + length(diag(A)) + length(diag(A)) + length(diag(A)))])
  }

  if (A.matrix == "full" & R.matrix == "symmetric") {

    A <- diag(c(init.par[1:m]))
    nr.off.diag <- upper.tri(A)
    l.upp.tri   <- length(nr.off.diag[nr.off.diag == TRUE])
    A[upper.tri(A)] <- init.par[(length(diag(A)) + 1):(length(diag(A)) + l.upp.tri)]
    A[lower.tri(A)] <- init.par[(length(diag(A)) + l.upp.tri + 1):(length(diag(A)) + l.upp.tri + l.upp.tri)]
    P <- eigen(A)$vectors
    D <- diag(eigen(A)$values)

    Chol <- diag(init.par[(length(diag(A)) + l.upp.tri + l.upp.tri + 1):(length(diag(A)) + l.upp.tri + l.upp.tri + length(diag(A)))])
    nr.off.diag.R <- upper.tri(Chol)
    l.upp.tri.R   <- length(nr.off.diag.R[nr.off.diag.R == TRUE])
    Chol[upper.tri(Chol)] <- init.par[(length(diag(A)) + l.upp.tri + l.upp.tri + length(diag(A)) + 1):(length(diag(A)) + l.upp.tri + l.upp.tri + length(diag(A)) + l.upp.tri.R)]

    optima <- c(init.par[(length(diag(A)) + l.upp.tri + l.upp.tri + length(diag(A)) + l.upp.tri.R + 1):(length(diag(A)) + l.upp.tri + l.upp.tri + length(diag(A)) + l.upp.tri.R + length(diag(A)))])
    anc    <- c(init.par[(length(diag(A)) + l.upp.tri + l.upp.tri + length(diag(A)) + l.upp.tri.R + length(diag(A)) + 1):(length(diag(A)) + l.upp.tri + l.upp.tri + length(diag(A)) + l.upp.tri.R + length(diag(A)) + length(diag(A)))])
  }

  if (A.matrix == "upper.tri" & R.matrix == "diag") {

    A <- diag(c(init.par[1:m]))
    nr.off.diag <- upper.tri(A)
    l.upp.tri   <- length(nr.off.diag[nr.off.diag == TRUE])
    A[upper.tri(A)] <- init.par[(length(diag(A)) + 1):(length(diag(A)) + l.upp.tri)]
    P <- eigen(A)$vectors
    D <- diag(eigen(A)$values)

    Chol   <- diag(init.par[(length(diag(A)) + l.upp.tri + 1):((length(diag(A)) + l.upp.tri + m))])
    optima <- c(init.par[(length(diag(A)) + l.upp.tri + m + 1):((length(diag(A)) + l.upp.tri + m + m))])
    anc    <- c(init.par[(length(diag(A)) + l.upp.tri + m + m + 1):(length(diag(A)) + l.upp.tri + m + m + m)])
  }

  if (A.matrix == "upper.tri" & R.matrix == "symmetric") {

    A <- diag(c(init.par[1:m]))
    nr.off.diag.A <- upper.tri(A)
    l.upp.tri.A   <- length(nr.off.diag.A[nr.off.diag.A == TRUE])
    A[upper.tri(A)] <- init.par[(length(diag(A)) + 1):(length(diag(A)) + l.upp.tri.A)]
    P <- eigen(A)$vectors
    D <- diag(eigen(A)$values)

    Chol <- diag(init.par[(length(diag(A)) + l.upp.tri.A + 1):(length(diag(A)) + l.upp.tri.A + m)])
    nr.off.diag.R <- upper.tri(Chol)
    l.upp.tri.R   <- length(nr.off.diag.R[nr.off.diag.R == TRUE])
    Chol[upper.tri(Chol)] <- init.par[(length(diag(A)) + l.upp.tri.A + m + 1):(length(diag(A)) + l.upp.tri.A + m + l.upp.tri.R)]

    optima <- c(init.par[(length(diag(A)) + l.upp.tri.A + m + l.upp.tri.R + 1):(length(diag(A)) + l.upp.tri.A + m + l.upp.tri.R + m)])
    anc    <- c(init.par[(length(diag(A)) + l.upp.tri.A + m + l.upp.tri.R + m + 1):(length(diag(A)) + l.upp.tri.A + m + l.upp.tri.R + m + m)])
  }

  if (A.matrix == "lower.tri" & R.matrix == "diag") {

    A <- diag(c(init.par[1:m]))
    nr.off.diag <- lower.tri(A)
    l.low.tri   <- length(nr.off.diag[nr.off.diag == TRUE])
    A[lower.tri(A)] <- init.par[(length(diag(A)) + 1):(length(diag(A)) + l.low.tri)]
    P <- eigen(A)$vectors
    D <- diag(eigen(A)$values)

    Chol   <- diag(init.par[(length(diag(A)) + l.low.tri + 1):((length(diag(A)) + l.low.tri + m))])
    optima <- c(init.par[(length(diag(A)) + l.low.tri + m + 1):((length(diag(A)) + l.low.tri + m + m))])
    anc    <- c(init.par[(length(diag(A)) + l.low.tri + m + m + 1):(length(diag(A)) + l.low.tri + m + m + m)])
  }

  if (A.matrix == "lower.tri" & R.matrix == "symmetric") {

    A <- diag(c(init.par[1:m]))
    nr.off.diag.A <- lower.tri(A)
    l.low.tri.A   <- length(nr.off.diag.A[nr.off.diag.A == TRUE])
    A[lower.tri(A)] <- init.par[(length(diag(A)) + 1):(length(diag(A)) + l.low.tri.A)]
    P <- eigen(A)$vectors
    D <- diag(eigen(A)$values)

    Chol <- diag(init.par[(length(diag(A)) + l.low.tri.A + 1):(length(diag(A)) + l.low.tri.A + m)])
    nr.off.diag.R <- upper.tri(Chol)
    l.upp.tri.R   <- length(nr.off.diag.R[nr.off.diag.R == TRUE])
    Chol[upper.tri(Chol)] <- init.par[(length(diag(A)) + l.low.tri.A + m + 1):(length(diag(A)) + l.low.tri.A + m + l.upp.tri.R)]

    optima <- c(init.par[(length(diag(A)) + l.low.tri.A + m + l.upp.tri.R + 1):(length(diag(A)) + l.low.tri.A + m + l.upp.tri.R + m)])
    anc    <- c(init.par[(length(diag(A)) + l.low.tri.A + m + l.upp.tri.R + m + 1):(length(diag(A)) + l.low.tri.A + m + l.upp.tri.R + m + m)])
  }

  if (A.matrix == "OUBM") {

    A <- diag(c(init.par[1:(m - 1)], 0.000001))
    nr.off.diag.A <- upper.tri(A)
    l.upp.tri.A   <- length(nr.off.diag.A[nr.off.diag.A == TRUE])
    A[upper.tri(A)] <- init.par[(m):(m - 1 + l.upp.tri.A)]
    P <- eigen(A)$vectors
    D <- diag(eigen(A)$values)

    Chol   <- diag(init.par[(m + l.upp.tri.A):(m + m + l.upp.tri.A - 1)])
    optima <- c(init.par[(m + m + l.upp.tri.A):(m + m + m + l.upp.tri.A - 1)])
    anc    <- c(init.par[(m + m + m + l.upp.tri.A):(m + m + m + m + l.upp.tri.A - 1)])
  }

### v.1.4 ###
# v.1.3  computed the time vector, the per-time-point exp(-A*t) array, M, and VV3 by rebuilding everything from `yy$tt` inside this 
# function on every call, using explicit i/j/k/l loops throughout. Rewritten below to use the pre-computed ta/tij/time_vec/se_vec arguments, 
# vectorised matrix expressions in place of the innermost loops, and to add tryCatch/finite-value guards (see end of function) that were 
#  not present in the release version. Also fixes a variable-shadowing bug in the previous version of the code where the loop variable `m` was reused inside 
#  the trait-count-consuming exponent loop
  
  ### -----------------------------------------------------------------------
  ### Optimization 1: pre-compute constants that are reused across time pairs
  ### -----------------------------------------------------------------------

  # P_inv: computed once here instead of inside the i,j loop (original called
  # solve(P) inside the loop, repeating the inversion n_t^2 times per call)
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
  ### Calculate expected trait values M (unified loop for all A.matrix types)
  ### -----------------------------------------------------------------------

  # In version 1.3, for diagonal A, the original used a fully vectorized formula across all time
  # points simultaneously; for non-diagonal A it used an explicit loop.  The
  # loop below is equivalent for all model types and also avoids building the
  # exp_eigenvalues array.

  anc_minus_theta <- anc - optima
  M <- matrix(NA, nrow = m, ncol = n_t)

  for (i in 1:n_t) {
    eAt_i  <- P %*% diag(exp(-d * time_vec[i])) %*% P_inv
    M[, i] <- eAt_i %*% anc_minus_theta + optima
  }
  M <- c(t(M))   # vectorize row-by-row to match ordering of y

  ### -----------------------------------------------------------------------
  ### Compute variance-covariance matrix VV3 (eq. 8 and 9 from Suppl. of Clavel et al. 2015)
  ###
  ### v.1.4 vs. v.1.3:
  ###   - No tmp.VV intermediate array; VV3 is filled directly.
  ###   - left.side computed by outer() instead of a k,l double loop.
  ###   - exp decay vector computed by a single vectorized exp() call
  ###     (also fixes variable-shadowing bug: original used 'm' as loop
  ###     variable, overwriting the trait-count variable).
  ###   - row/col block offsets pre-computed outside the inner loop.
  ### -----------------------------------------------------------------------

  VV3 <- matrix(0, nrow = n_t * m, ncol = n_t * m)

  # Block-row and block-column starting offsets for each trait in VV3.
  # VV3[(trait_r - 1)*n_t + time_a, (trait_c - 1)*n_t + time_b] holds the
  # covariance of trait c at time a with trait r at time b (see assembly logic).
  row_base <- (0:(m - 1)) * n_t   # length-m vector of row offsets
  col_base <- row_base             # same for columns (symmetric block structure)

  for (a in 1:n_t) {
    for (b in 1:n_t) {

      # Vectorized left.side: replaces the original double loop over k and l.
      left.side  <- (1 - exp(-sum_d * ta[a, b])) / sum_d   # m x m matrix
      left.right <- left.side * right.side                  # Hadamard product
      integ      <- P %*% left.right %*% t(P)

      # Vectorized exp decay: replaces the original inner loop over traits.
      # Variable name 'exp_diag_b' avoids shadowing the outer 'm'.
      exp_diag_b <- exp(-d * tij[a, b])                     # length-m vector
      exp_mat_T  <- t(P %*% diag(exp_diag_b) %*% P_inv)    # transpose of exp(-A*tij)

      tmp <- integ %*% exp_mat_T   # m x m block for time pair (a, b)

      # Direct placement into VV3.
      # VV3[row_base + a, col_base + b] is the m x m submatrix corresponding to
      # time point a (rows) crossed with time point b (columns) across all trait
      # pairs.  Setting it to t(tmp) replicates the block assembly from the
      # original tmp.VV / List approach.
      row_idx <- row_base + a
      col_idx <- col_base + b
      VV3[row_idx, col_idx] <- t(tmp)
    }
  }

  ### Add pre-computed sampling error to the diagonal ###
  diag(VV3) <- diag(VV3) + se_vec
  VV3 <- (VV3 + t(VV3)) / 2   # guard against floating-point asymmetry

  if (any(!is.finite(M)) || any(!is.finite(VV3))) return(-1e20)

  S <- tryCatch(
    mvtnorm::dmvnorm(t(y), mean = M, sigma = VV3, log = TRUE),
    error = function(e) -1e20
  )
  if (!is.finite(S)) return(-1e20)
  S
}
