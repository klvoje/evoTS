#' @title Log-likelihoods for evolutionary models
#'
#' @description Returns log-likelihood for a model with an Unbiased Random Walk in the first segment and an Unbiased Random Walk in the second segment.
#'
#' @param p parameters of the model to be optimized
#'
#' @param y a paleoTS object
#'
#' @param gg numeric vector indicating membership of each sample in a segment
#'
#' @details In general, users will not access these functions directly, but instead use the optimization functions, which use these functions to find the best-supported parameter values.
#'
#'@return The log-likelihood of the parameter estimates, given the data.
#'
#'@author Kjetil Lysne Voje


logL.joint.URW.URW<-function (p, y, gg)
{
  #These parameters need to mach the order of p0
  anc<-p[1]
  vs_1 <- p[2]
  vs_2 <- p[3]

  n <- length(y$mm)
  #creates a vector that repeates the ancestral value (first sample mean in sequence)
  M <- rep(anc, n)
  M <- unname(M)

  #define a vector based on length of first and second RW in the time series
  rw_1.seg <- which(gg == 1)
  rw_2.seg <- which(gg == 2)

  #Extract the time vector for the random walks:
  tt.rw_1 <- y$tt[rw_1.seg]
  #tt.rw_2 <- y$tt[rw_2.seg] - y$tt[rw_2.seg[1] - 1]
  tt.rw_2 <- y$tt[rw_2.seg] - y$tt[rw_2.seg[1]]

  #tt.rw_2 <- y$tt[rw_2.seg]

  #Multiply variance parameter (vstep) from RW with outer product of time*time
  VVrw_1 <- vs_1 * outer(tt.rw_1, tt.rw_1, FUN = pmin)
  VVrw_2 <- vs_2 * outer(tt.rw_2, tt.rw_2, FUN = pmin)

  #Create empty n*n matrix
  VVtot <- array(0, dim = c(n, n))

  #Create final variance matrix by combining the variance matrices from stasis and RW
  VVtot[rw_1.seg, rw_1.seg] <- VVrw_1
  VVtot[rw_2.seg, rw_2.seg] <- VVrw_2

  #Add population variance to the diagonal of the variance matrix
  diag(VVtot) <- diag(VVtot) + y$vv/y$nn

  #
  S <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VVtot, log = TRUE)
  return(S)
}
