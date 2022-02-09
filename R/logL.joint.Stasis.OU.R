#' @title Log-likelihoods for evolutionary models
#'
#' @description Returns log-likelihood for a model with stasis in the first segment and an Ornstein-Uhlenbeck process in the second segment. .
#'
#' @param p parameters of the model to be optimized
#'
#' @param y a paleoTS object
#'
#' @param gg numeric vector indicating membership of each sample in a segment
#'
#' @details In general, users will not be access these functions directly, but instead use the optimization functions, which use these functions to find the best-supported parameter values.
#'
#'@return The log-likelihood of the parameter estimates, given the data.
#'
#'@author Kjetil Lysne Voje


logL.joint.Stasis.OU<-function (p, y, gg)
{
  #These parameters need to match the order of p0
  theta_st <- p[1]
  omega <- p[2]
  vstep <- p[3]
  theta_OU <- p[4]
  alpha <- p[5]

  n <- length(y$mm)

  #define a vector based on length of first and second RW in the time series
  st.seg <- which(gg == 1)
  OU.seg <- which(gg == 2)

  #Extract the time vector for the OU:
  tt.OU <- y$tt[OU.seg] - y$tt[OU.seg[1] - 1]

  #create a vector for the expected means for the whole time series
  M <- c(rep(theta_st, length(st.seg)), ou.M(theta_st, theta_OU, alpha, tt.OU))
  M <- unname(M)

  #Create variance-covariance matrices for the three models
  VV.st <- diag(omega, nrow = length(st.seg))

  ff <- function(a, b) abs(a - b)
  VV.OU <- outer(tt.OU, tt.OU, FUN = ff)
  VV.OU <- exp(-alpha * VV.OU)
  VVd.OU <- ou.V(vstep, alpha, tt.OU)
  VV2.OU <- outer(VVd.OU, VVd.OU, pmin)
  VV.OU <- VV.OU * VV2.OU

  #Create empty n*n matrix
  VVtot <- array(0, dim = c(n, n))

  #Create final variance matrix by combining the variance matrices from stasis and RW
  VVtot[st.seg, st.seg] <- VV.st
  VVtot[OU.seg, OU.seg] <- VV.OU

  #Add population variance to the diagonal of the variance matrix
  diag(VVtot) <- diag(VVtot) + y$vv/y$nn

  #
  S <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VVtot, log = TRUE)
  return(S)
}
