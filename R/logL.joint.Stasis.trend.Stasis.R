#' @title Log-likelihoods for evolutionary models
#'
#' @description Returns log-likelihood for a model with stasis in the first segment, directional trend (GRW), and stasis in the third segment.
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


logL.joint.Stasis.trend.Stasis<-function (p, y, gg)
{
  #These parameters need to match the order of p0
  theta_1 <- p[1]
  omega_1 <- p[2]
  mstep <- p[3]
  vstep <- p[4]
  theta_2 <- p[5]
  omega_2 <- p[6]
  n <- length(y$mm)

  #define a vector based on length of first and second RW in the time series
  st_1.seg <- which(gg == 1)
  st_2.seg <- which(gg == 3)
  dt.seg <- which(gg == 2)

  #Extract the time vector for the trend:
  tt.trend <- y$tt[dt.seg] - y$tt[dt.seg[1]]

  #create a vector for the expected means for the whole time series
  M <- c(rep(theta_1, length(st_1.seg)), theta_1 + mstep * tt.trend, rep(theta_2, length(st_2.seg)))
  M <- unname(M)

  #Create variance-covariance matrices for the three models
  VVst_1 <- diag(omega_1, nrow = length(st_1.seg))
  VVtrend <- vstep * outer(tt.trend, tt.trend, FUN = pmin)
  VVst_2 <- diag(omega_2, nrow = length(st_2.seg))

  #Create empty n*n matrix
  VVtot <- array(0, dim = c(n, n))

  #Create final variance matrix by combining the variance matrices from stasis and RW
  VVtot[st_1.seg, st_1.seg] <- VVst_1
  VVtot[dt.seg, dt.seg] <- VVtrend
  VVtot[st_2.seg, st_2.seg] <- VVst_2

  #Add population variance to the diagonal of the variance matrix
  diag(VVtot) <- diag(VVtot) + y$vv/y$nn

  #
  S <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VVtot, log = TRUE)
  return(S)
}
