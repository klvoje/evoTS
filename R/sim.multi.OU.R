#' @title Simulate multivariate Ornstein-Uhlenbeck evolutionary sequence data sets
#'
#' @description Function to simulate a multivariate Ornstein-Uhlenbeck evolutionary sequence data set.
#'
#' @param ns number of samples in time-series
#'
#' @param anc the ancestral trait values
#'
#' @param optima the optimal trait values
#'
#' @param A the pull matrix.
#'
#' @param R the drift matrix
#'
#' @param vp within-population trait variance
#'
#' @param nn 	vector of the number of individuals in each sample (identical sample sizes for all time-series is assumed)
#'
#' @param tt 	vector of sample ages, increases from oldest to youngest
#'
#'@return A multivariate evolutionary sequence (time-series) data set.
#'
#'@note The Ornstein Uhlenbeck model is reduced to an Unbiased Random Walk when the alpha parameter is zero. It is therefore possible to let a trait evolve as an Unbiased Random Walk by setting the diagonal element for that trait to a value close to zero (e.g. 1e-07). Elements in the diagonal of A cannot be exactly zero as this will result in a singular variance-covariance matrix.
#'
#'@author Kjetil Lysne Voje
#'
#'@export
#'
#'@examples
#'##Define the A and R matrices
#'
#'A_matrix<-matrix(c(4,-2,0,3), nrow=2, byrow = TRUE)
#'R_matrix<-matrix(c(4,0.2,0.2,4), nrow=2, byrow = TRUE)
#'
#'## Generate an evoTS object by simulating a multivariate dataset
#'data_set<-sim.multi.OU(40, optima = c(1.5,2),A=A_matrix , R = R_matrix)
#'
#'## plot the data
#'plotevoTS.multivariate(data_set)
#'
sim.multi.OU<-function(ns = 30, anc = c(0,0), optima = c(3, 2),
                       A= matrix(c(7,0,0,6), nrow=2, byrow = TRUE), R = matrix(c(0.05,0,0,0.05), nrow=2, byrow = TRUE),
                       vp = 0.1, nn = rep(30, ns), tt = 0:(ns - 1)){
  m<-ncol(A)

### v.1.4 ###
# An almost tottal rewrite of the codebody.
#
# 1) Two bugs are fixed:
#    a) V.1.3's observed-value sampling loop was
#       `for (i in 2:ns) { x <- MASS::mvrnorm(nn[j], mu = MM[j,i], Sigma = sqrt(vp)); ... }`
#       (see below) — it started at i = 2, so mm[,1]/vv[,1] were never
#       assigned and stayed NA for the first time point. It also indexed
#       sample size as nn[j] (a trait index) instead of nn[i] (a time index).
#       Both are fixed below (loop now runs 1:ns, uses nn[i], and uses plain
#       rnorm() instead of the mvrnorm()-for-a-scalar workaround the release
#       version used).
#    b) The release version's expected-trajectory calculation
#       (`traits[,i]<-((P%*%diag(c(exp(-diag(D)[1]*time[1,i]),exp(-diag(D)[2]*time[1,i])))%*%solve(P))%*%anc) + ...`)
#       was hardcoded for exactly 2 traits. Generalized below to arbitrary m
#       via exp_At <- P %*% diag(exp(-d * time[i])) %*% P.inv.
#
# 2) Performance: P.inv, right.side, and d_sum are now pre-computed once
#    instead of being recomputed inside the i,j double loop; the inner k,l
#    loop that built left.side element-by-element is replaced by a single
#    vectorised matrix expression.
#

  A.matrix <- A
  P        <- eigen(A.matrix)$vectors
  D        <- diag(eigen(A.matrix)$values)
  d        <- diag(D)      # eigenvalues as vector
  P.inv    <- solve(P)     # pre-compute once; used repeatedly inside loops
  Chol     <- chol(R)

  # Normalised time vector
  time <- tt / max(tt)

  tij <- outer(time, time, FUN = function(a, b) abs(a - b))
  ta  <- outer(time, time, FUN = pmin)

  # Pre-compute right side of the VCV integral (constant across all time pairs)
  right.side <- P.inv %*% (Chol %*% t(Chol)) %*% t(P.inv)

  # Pre-compute eigenvalue sum matrix for vectorised left.side computation
  d_sum <- outer(d, d, "+")

  # Compute expected trait values at each time point (generalised to any m)
  traits <- matrix(NA, nrow = m, ncol = ns)
  for (i in 1:ns) {
    exp_At     <- P %*% diag(exp(-d * time[i])) %*% P.inv
    traits[,i] <- exp_At %*% anc + (diag(m) - exp_At) %*% optima
  }

  # Build the VCV array: tmp.VV[i,j,] is the vectorised m x m block for time pair (i,j)
  tmp.VV <- array(data = NA, dim = c(ns, ns, m * m))

  for (i in 1:ns) {
    for (j in 1:ns) {
      # Vectorised replacement of the inner k,l loops
      left.side  <- (1 - exp(-d_sum * ta[i,j])) / d_sum
      left.right <- left.side * right.side
      integ      <- P %*% left.right %*% t(P)

      exp_eigenvalues <- exp(-d * tij[i,j])
      tmp             <- integ %*% t(P %*% diag(exp_eigenvalues) %*% P.inv)
      tmp.VV[i, j, ]  <- as.vector(tmp)
    }
  }

  # Assemble the full block VCV matrix
  VV3 <- matrix(0, ncol = ns * m, nrow = ns * m)
  List <- vector("list", m * m)
  for (k in 1:(m * m)) List[[k]] <- tmp.VV[,,k]

  from.boundary <- seq(1, ns * m, ns)
  to.boundary   <- from.boundary + ns - 1
  from          <- seq(1, m * m, m)
  to            <- from + m - 1

  for (i in 1:m) {
    VV3[from.boundary[i]:to.boundary[i], ] <- do.call(cbind, List[from[i]:to[i]])
  }

  # Draw population means from multivariate normal
  traits_x <- c(t(traits))
  MM_vec   <- MASS::mvrnorm(n = 1, mu = traits_x, Sigma = VV3)
  MM       <- matrix(MM_vec, m, ns, byrow = TRUE)

  # Sample observed means and within-sample variances
  mm <- matrix(nrow = m, ncol = ns)
  vv <- matrix(nrow = m, ncol = ns)

  for (i in 1:ns) {
    for (j in 1:m) {
      x       <- rnorm(nn[i], mean = MM[j,i], sd = sqrt(vp))
      mm[j,i] <- mean(x)
      vv[j,i] <- var(x)
    }
  }

  List2<-list()
  for (i in 1:m){
    List2[[i]]<-paleoTS::as.paleoTS(mm = mm[i,], vv = vv[i,], nn = nn, tt = time, MM = MM[i,], label = "Created by sim.multi.OU()", reset.time = FALSE)
  }

  if (m==2) yy<-make.multivar.evoTS(List2[[1]], List2[[2]])
  if (m==3) yy<-make.multivar.evoTS(List2[[1]], List2[[2]], List2[[3]])
  if (m==4) yy<-make.multivar.evoTS(List2[[1]], List2[[2]], List2[[3]], List2[[4]])
  if (m==5) yy<-make.multivar.evoTS(List2[[1]], List2[[2]], List2[[3]], List2[[4]], List2[[5]])

  return(yy)
}
