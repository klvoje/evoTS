#' @title Log-likelihoods for evolutionary models
#'
#' @description Returns log-likelihood for a model with an Ornstein-Uhlenbeck process in the first segment, an Ornstein-Uhlenbeck process in the second segment, and an Ornstein-Uhlenbeck process in the third segment.
#'
#' @details In general, users will not access these functions directly, but instead use the optimization functions, which use these functions to find the best-supported parameter values.
#'
#' @param p initial (starting) parameter values
#'
#' @param y a paleoTS object
#'
#' @param gg numeric vector indicating membership of each sample in a segment
#'
#'@return The log-likelihood of the parameter estimates, given the data.
#'
#'@author Kjetil Lysne Voje

logL.joint.OU.OU.OU<-function (p, y, gg)
{
  #These parameters need to match the order of p0
  anc<-      p[1]
  vstep_1 <- p[2]
  theta_1 <- p[3]
  alpha_1 <- p[4]
  vstep_2 <- p[5]
  theta_2 <- p[6]
  alpha_2 <- p[7]
  vstep_3 <- p[8]
  theta_3 <- p[9]
  alpha_3 <- p[10]

  n <- length(y$mm)

  ou.M<-function (anc, theta, aa, tt) { theta * (1 - exp(-aa * tt)) + anc * exp(-aa * tt)}
  ou.V<-function (vs, aa, tt) { (vs/(2 * aa)) * (1 - exp(-2 * aa * tt))}

  #define a vector based on length of first and second RW in the time series
  OU.seg_1 <- which(gg == 1)
  OU.seg_2 <- which(gg == 2)
  OU.seg_3 <- which(gg == 3)


  #Extract the time vector for the OU:
  tt.OU_1 <- y$tt[OU.seg_1]
  tt.OU_2 <- y$tt[OU.seg_2] - y$tt[OU.seg_2[1]]
  tt.OU_3 <- y$tt[OU.seg_3] - y$tt[OU.seg_3[1]]


  #create a vector for the expected means for the whole time series
  M <- c(ou.M(anc, theta_1, alpha_1, tt.OU_1), ou.M(theta_1, theta_2, alpha_2, tt.OU_2),ou.M(theta_2, theta_3, alpha_3, tt.OU_3))
  M <- unname(M)

  #Create variance-covariance matrices for the three models
  ff <- function(a, b) abs(a - b)

  VV.OU_1 <- outer(tt.OU_1, tt.OU_1, FUN = ff)
  VV.OU_1 <- exp(-alpha_1 * VV.OU_1)
  VVd.OU_1 <- ou.V(vstep_1, alpha_1, tt.OU_1)
  VV2.OU_1 <- outer(VVd.OU_1, VVd.OU_1, pmin)
  VV.OU_1 <- VV.OU_1 * VV2.OU_1

  VV.OU_2 <- outer(tt.OU_2, tt.OU_2, FUN = ff)
  VV.OU_2 <- exp(-alpha_2 * VV.OU_2)
  VVd.OU_2 <- ou.V(vstep_2, alpha_2, tt.OU_2)
  VV2.OU_2 <- outer(VVd.OU_2, VVd.OU_2, pmin)
  VV.OU_2 <- VV.OU_2 * VV2.OU_2

  VV.OU_3 <- outer(tt.OU_3, tt.OU_3, FUN = ff)
  VV.OU_3 <- exp(-alpha_3 * VV.OU_3)
  VVd.OU_3 <- ou.V(vstep_3, alpha_3, tt.OU_3)
  VV2.OU_3 <- outer(VVd.OU_3, VVd.OU_3, pmin)
  VV.OU_3 <- VV.OU_3 * VV2.OU_3

  #Create empty n*n matrix
  VVtot <- array(0, dim = c(n, n))

  #Create final variance matrix by combining the variance matrices from stasis and RW
  VVtot[OU.seg_1, OU.seg_1] <- VV.OU_1
  VVtot[OU.seg_2, OU.seg_2] <- VV.OU_2
  VVtot[OU.seg_3, OU.seg_3] <- VV.OU_3


  #Add population variance to the diagonal of the variance matrix
  diag(VVtot) <- diag(VVtot) + y$vv/y$nn

  #
  S <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VVtot, log = TRUE)
  return(S)
}
