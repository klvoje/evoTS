#' @title Optimization and log-likelihoods for pairs of models.
#'
#' @description A collections of functions that serves the function fit.mode.shift. See fit.mode.shift for info.
#'
#' @param y a paleoTS object.
#'
#' @param gg numeric vector indicating membership of each sample in segments
#'
#' @param cl control list to be passed to optim
#'
#' @param pool logical indicating whether to pool variances across samples
#'
#' @param meth optimization method, passed to function optim. Default is "L-BFGS-B".
#'
#' @param hess logical, indicating whether to calculate standard errors from the Hessian matrix.
#'
#' @details In general, users will not be access these functions directly, but instead use the wrapper function, which use these functions to find the best-supported parameter values.
#'
#' @note This function is not likely to be called directly by the user.
#'
#' @author Kjetil Lysne Voje
#'

opt.joint.URW.Stasis<-function (y, gg, cl = list(fnscale = -1), pool = TRUE, meth = "L-BFGS-B", hess = FALSE)
{
  if (pool) {
    y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  }
  small <- 1e-08
  p0st <- paleoTS::mle.Stasis(paleoTS::sub.paleoTS(y, ok = gg == 2))
  names(p0st) <- c("theta", "omega")
  p0urw <- paleoTS::mle.URW(paleoTS::sub.paleoTS(y, ok = gg == 1))
  names(p0urw) <-"vstep"
  p0anc<-y$mm[1]
  names(p0anc)<-"anc"

  K <- 5

  #Checking if the variance parameter in the RW and Stasis is below critical value, and if TRUE, multiply it with 100
  if (p0st["omega"] <= small) p0st["omega"] <- 100 * small
  if (p0urw["vstep"] <= small) p0urw["vstep"] <- 100 * small

  #make vector out of estimated model parameters
  p0 <- c(p0anc, p0urw, p0st)

  #define lower bounds for the L-BFGS-B method
  # three parameters in the unbiased ramdom walk and 4 parameters for the biased random walk
  ll<- c(NA, small, NA, small)

  # If user wants to use a differen method than L-BFGS-B, then the lower boun is defined as -Inf
  if (meth != "L-BFGS-B")
    ll <- -Inf

  #Here the multivariate parameter estimation routine start:
  w <- try(optim(p0, fn = logL.joint.URW.Stasis, gg = gg,
                 method = meth, lower = ll, control = cl, hessian = hess,
                 y = y), silent = TRUE)
  if (inherits (w, "try-error")) {
    cl <- list(fnscale = -1, parscale = c(1, 10, 100))
    w <- try(optim(p0, fn = logL.joint.URW.Stasis, gg = gg,
                   method = meth, lower = ll, control = cl, hessian = hess,
                   y = y), silent = TRUE)
  }

  if (inherits (w, "try-error")) {
    wc <- paleoTS::as.paleoTSfit(logL = NA, parameters = NA, modelName = paste("URW",
                                                                      "Stasis", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                        se = NULL)
    return(wc)
  }
  if (hess)
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- paleoTS::as.paleoTSfit(logL = w$value, parameters = w$par, modelName = paste("URW",
                                                                            "Stasis", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                      se = w$se)
  return(wc)
}

logL.joint.URW.Stasis<-function (p, y, gg)
{
  anc   <- p[1]
  vstep <- p[2]
  theta <- p[3]
  omega <- p[4]

  n <- length(y$mm)
  M <- rep(anc, n)

  #define a vector based on length of first and second RW in the time series
  st.seg <- which(gg == 2)
  urw.seg <- which(gg == 1)
  tt.urw <- y$tt[urw.seg] - y$tt[urw.seg[1]]
  #tt.urw <- y$tt[urw.seg]

  #Create variance-covariance matrices for the three models
  VV.st <- diag(omega, nrow = length(st.seg))
  VV.urw <- vstep * outer(tt.urw, tt.urw, FUN = pmin)

  #Create empty n*n matrix
  VVtot <- array(0, dim = c(n, n))

  #Create final variance matrix by combining the variance matrices from stasis and RW
  VVtot[st.seg, st.seg] <- VV.st
  VVtot[urw.seg, urw.seg] <- VV.urw

  #Add population variance to the diagonal of the variance matrix
  diag(VVtot) <- diag(VVtot) + y$vv/y$nn

  #
  S <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VVtot, log = TRUE)
  return(S)
}

########

opt.joint.GRW.Stasis<-function (y, gg, cl = list(fnscale = -1), pool = TRUE, meth = "L-BFGS-B", hess = FALSE)
{
  if (pool) {
    y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  }
  small <- 1e-08
  p0st <- paleoTS::mle.Stasis(paleoTS::sub.paleoTS(y, ok = gg == 2))
  p0dt<- paleoTS::mle.GRW(paleoTS::sub.paleoTS(y, ok = gg == 1))
  p0anc<-y$mm[1]
  names(p0anc)<-"anc"

  K <- 6

  #Checking if te variance parameter in the RW and Stasis is below critical value, and if TRUE, multiply it with 100
  if (p0st["omega"] <= small) p0st["omega"] <- 100 * small
  if (p0dt["vstep"] <= small) p0dt["vstep"] <- 100 * small

  #make vector out of estimated model parameters
  p0 <- c(p0anc, p0dt, p0st)

  #define lower bounds for the L-BFGS-B method
  # three parameters in the unbiased ramdom walk and 4 parameters for the biased random walk
  ll<- c(NA, NA, small, NA, small)

  # If user wants to use a differen method than L-BFGS-B, then the lower boun is defined as -Inf
  if (meth != "L-BFGS-B")
    ll <- -Inf

  #Here the multivariate parameter estimation routine start:
  w <- try(optim(p0, fn = logL.joint.GRW.Stasis, gg = gg,
                 method = meth, lower = ll, control = cl, hessian = hess,
                 y = y), silent = TRUE)
  if (inherits (w, "try-error")) {
    cl <- list(fnscale = -1, parscale = c(1, 10, 100))
    w <- try(optim(p0, fn = logL.joint.GRW.Stasis, gg = gg,
                   method = meth, lower = ll, control = cl, hessian = hess,
                   y = y), silent = TRUE)
  }

  if (inherits (w, "try-error")) {
    wc <- paleoTS::as.paleoTSfit(logL = NA, parameters = NA, modelName = paste("GRW",
                                                                      "Stasis", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                        se = NULL)
    return(wc)
  }
  if (hess)
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- paleoTS::as.paleoTSfit(logL = w$value, parameters = w$par, modelName = paste("GRW",
                                                                            "Stasis", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                      se = w$se)
  return(wc)
}

logL.joint.GRW.Stasis<-function (p, y, gg)
{
  #These parameters need to match the order of p0
   anc  <- p[1]
  mstep <- p[2]
  vstep <- p[3]
  theta <- p[4]
  omega <- p[5]

  n <- length(y$mm)

  #define a vector based on length of first and second RW in the time series
  st.seg <- which(gg == 2)
  dt.seg <- which(gg == 1)

  #Extract the time vector for the trend:
  tt.trend <- y$tt[dt.seg]


  #create a vector for the expected means for the whole time series
  M <- c(anc + mstep * tt.trend, rep(theta, length(st.seg)))
  M <- unname(M)

  #Create variance-covariance matrices for the three models
  VVst <- diag(omega, nrow = length(st.seg))
  VVtrend <- vstep * outer(tt.trend, tt.trend, FUN = pmin)

  #Create empty n*n matrix
  VVtot <- array(0, dim = c(n, n))

  #Create final variance matrix by combining the variance matrices from stasis and RW
  VVtot[st.seg, st.seg] <- VVst
  VVtot[dt.seg, dt.seg] <- VVtrend

  #Add population variance to the diagonal of the variance matrix
  diag(VVtot) <- diag(VVtot) + y$vv/y$nn

  S <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VVtot, log = TRUE)
  return(S)
}

######

opt.joint.GRW.URW<-function (y, gg, cl = list(fnscale = -1), pool = TRUE, meth = "L-BFGS-B", hess = FALSE)
{
  if (pool) {
    y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  }
  small <- 1e-08
  p0rw <- paleoTS::mle.URW(paleoTS::sub.paleoTS(y, ok = gg == 2))
  p0dt <- paleoTS::mle.GRW(paleoTS::sub.paleoTS(y, ok = gg == 1))

  K <- 5

  #Checking if te variance parameter in the RW and Stasis is below critical value, and if TRUE, multiply it with 100
  if (p0rw["vstep"] <= small)
    p0rw["vstep"] <- 100 * small
  if (p0dt["vstep"] <= small)
    p0dt["vstep"] <- 100 * small

  # Define ancestral state
  p0anc <- y$mm[1]
  names(p0anc) <- "anc"
  names(p0rw)<-"vstep.urw"
  names(p0dt)<-c("mstep", "vstep.grw")

  #make vector out of estimated model parameters
  p0 <- c(p0anc, p0dt, p0rw)

  #define lower bounds for the L-BFGS-B method
  # three parameters in the unbiased ramdom walk and 4 parameters for the biased random walk
  ll<- c(NA, NA, small, small)

  # If user wants to use a differen method than L-BFGS-B, then the lower boun is defined as -Inf
  if (meth != "L-BFGS-B")
    ll <- -Inf

  #Here the multivariate parameter estimation routine start:
  w <- try(optim(p0, fn = logL.joint.GRW.URW, gg = gg,
                 method = meth, lower = ll, control = cl, hessian = hess,
                 y = y), silent = TRUE)
  if (inherits (w, "try-error")) {
    cl <- list(fnscale = -1, parscale = c(1, 10, 100))
    w <- try(optim(p0, fn = logL.joint.GRW.URW, gg = gg,
                   method = meth, lower = ll, control = cl, hessian = hess,
                   y = y), silent = TRUE)
  }

  if (inherits (w, "try-error")) {
    wc <- paleoTS::as.paleoTSfit(logL = NA, parameters = NA, modelName = paste("GRW",
                                                                      "URW", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                        se = NULL)
    return(wc)
  }
  if (hess)
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- paleoTS::as.paleoTSfit(logL = w$value, parameters = w$par, modelName = paste("GRW",
                                                                            "URW", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                      se = w$se)
  return(wc)
}

logL.joint.GRW.URW<-function (p, y, gg)
{
  #These parameters need to mach the order of p0
  anc<-p[1]
  mstep <- p[2]
  vstep.grw <- p[3]
  vstep.urw <- p[4]

  n <- length(y$mm)

  #define a vector based on length of first and second RW in the time series
  rw.seg <- which(gg == 2)
  dt.seg <- which(gg == 1)

  #Extract the time vector for the random walks:
  tt.rw <- y$tt[rw.seg] - y$tt[rw.seg[1]]
  #tt.rw <- y$tt[rw.seg]
  tt.dt <- y$tt[dt.seg]

  anc.2<-tail((anc + mstep * tt.dt),1)
  #creates a vector that repeates the ancestral value (first sample mean in sequence)
  M <- c((anc + mstep * tt.dt), rep(anc.2, length(tt.rw)))
  M <- unname(M)

  #Multiply variance parameter (vstep) from RW with outer product of time*time
  VVrw <- vstep.urw * outer(tt.rw, tt.rw, FUN = pmin)
  VVdt <- vstep.grw * outer(tt.dt, tt.dt, FUN = pmin)

  #Create empty n*n matrix
  VVtot <- array(0, dim = c(n, n))

  #Create final variance matrix by combining the variance matrices from stasis and RW
  VVtot[rw.seg, rw.seg] <- VVrw
  VVtot[dt.seg, dt.seg] <- VVdt

  #Add population variance to the diagonal of the variance matrix
  diag(VVtot) <- diag(VVtot) + y$vv/y$nn

  #
  S <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VVtot, log = TRUE)
  return(S)
}

#########

opt.joint.OU.Stasis<-function (y, gg, cl = list(fnscale = -1), pool = TRUE, meth = "L-BFGS-B", hess = FALSE)
{
  if (pool) {
    y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  }
  small <- 1e-08
  p0st <- paleoTS::mle.Stasis(paleoTS::sub.paleoTS(y, ok = gg == 2))[2]
  names(p0st) <- ("omega")
  p0OU<- paleoTS::mle.GRW(paleoTS::sub.paleoTS(y, ok = gg == 1))

  K <- 6

  #Checking if te variance parameter in the RW and Stasis is below critical value, and if TRUE, multiply it with 100
  if (p0st["omega"] <= small)
    p0st["omega"] <- 100 * small
  if (p0OU["vstep"] <= small)
    p0OU["vstep"] <- 100 * small

  #prepare initial guesses of OU parameters
  halft <- (paleoTS::sub.paleoTS(y, ok = gg == 1)$tt[length(paleoTS::sub.paleoTS(y, ok = gg == 1)$tt)] - paleoTS::sub.paleoTS(y, ok = gg == 1)$tt[1])/4
  p0 <- c(p0OU["vstep"]/10, paleoTS::sub.paleoTS(y, ok = gg == 1)$mm[length(paleoTS::sub.paleoTS(y, ok = gg == 1)$mm)], log(2)/halft)
  names(p0) <- c("vstep", "theta_OU", "alpha")

  p0anc<-y$mm[1]
  names(p0anc)<-"anc"
  #make vector out of estimated model parameters
  p0 <- c(p0anc, p0, p0st )

  #define lower bounds for the L-BFGS-B method
  # three parameters in the unbiased ramdom walk and 4 parameters for the biased random walk
  ll<- c(NA, small, NA, small, small)

  # If user wants to use a differen method than L-BFGS-B, then the lower boun is defined as -Inf
  if (meth != "L-BFGS-B")
    ll <- -Inf

  #Here the multivariate parameter estimation routine start:
  w <- try(optim(p0, fn = logL.joint.OU.Stasis, gg = gg,
                 method = meth, lower = ll, control = cl, hessian = hess,
                 y = y), silent = TRUE)
  if (inherits (w, "try-error")) {
    cl <- list(fnscale = -1, parscale = c(1, 10, 100))
    w <- try(optim(p0, fn = logL.joint.OU.Stasis, gg = gg,
                   method = meth, lower = ll, control = cl, hessian = hess,
                   y = y), silent = TRUE)
  }

  if (inherits (w, "try-error")) {
    wc <- paleoTS::as.paleoTSfit(logL = NA, parameters = NA, modelName = paste("OU",
                                                                      "Stasis", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                        se = NULL)
    return(wc)
  }
  if (hess)
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- paleoTS::as.paleoTSfit(logL = w$value, parameters = w$par, modelName = paste("OU",
                                                                            "Stasis", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                      se = w$se)
  return(wc)
}


logL.joint.OU.Stasis<-function (p, y, gg)
{
  #These parameters need to match the order of p0

  anc<-p[1]
  vstep <- p[2]
  theta_OU <- p[3]
  alpha <- p[4]
  omega <- p[5]

  n <- length(y$mm)

  #define a vector based on length of first and second RW in the time series
  st.seg <- which(gg == 2)
  OU.seg <- which(gg == 1)

  #Extract the time vector for the OU:
  tt.OU <- y$tt[OU.seg]

  #create a vector for the expected means for the whole time series
  M <- c(ou.M(anc, theta_OU, alpha, tt.OU), rep(theta_OU, length(st.seg)))
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

########

opt.joint.OU.URW<-function (y, gg, cl = list(fnscale = -1), pool = TRUE, meth = "L-BFGS-B", hess = FALSE)
{
  if (pool) {
    y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  }
  small <- 1e-08
  p0rw <- paleoTS::mle.URW(paleoTS::sub.paleoTS(y, ok = gg == 2))
  p0anc<-y$mm[1]
  names(p0anc) <- "anc"
  names(p0rw)<-"vstep.urw"
  halft <- (paleoTS::sub.paleoTS(y, ok = gg == 1)$tt[length(paleoTS::sub.paleoTS(y, ok = gg == 1)$tt)] - paleoTS::sub.paleoTS(y, ok = gg == 1)$tt[1])/4
  p0OU <- c(paleoTS::mle.GRW(paleoTS::sub.paleoTS(y, ok = gg == 1))[2]/10, paleoTS::sub.paleoTS(y, ok = gg == 1)$mm[length(paleoTS::sub.paleoTS(y, ok = gg == 1)$mm)], log(2)/halft)
  names(p0OU) <- c("vstep.OU", "theta_OU", "alpha")

  K <- 6

  #Checking if te variance parameter in the RW and Stasis is below critical value, and if TRUE, multiply it with 100
  if (p0rw["vstep.urw"] <= small) p0rw["vstep.urw"] <- 100 * small
  if (p0OU["vstep.OU"] <= small) p0OU["vstep.OU"] <- 100 * small

  #make vector out of estimated model parameters
  p0 <- c(p0anc, p0OU, p0rw)

  #define lower bounds for the L-BFGS-B method
  # three parameters in the unbiased ramdom walk and 4 parameters for the biased random walk
  ll<- c(NA, small, NA, small, small)

  # If user wants to use a differen method than L-BFGS-B, then the lower boun is defined as -Inf
  if (meth != "L-BFGS-B")
    ll <- -Inf

  #Here the multivariate parameter estimation routine start:
  w <- try(optim(p0, fn = logL.joint.OU.URW, gg = gg,
                 method = meth, lower = ll, control = cl, hessian = hess,
                 y = y), silent = TRUE)
  if (inherits (w, "try-error")) {
    cl <- list(fnscale = -1, parscale = c(1, 10, 100))
    w <- try(optim(p0, fn = logL.joint.OU.URW, gg = gg,
                   method = meth, lower = ll, control = cl, hessian = hess,
                   y = y), silent = TRUE)
  }

  if (inherits (w, "try-error")) {
    wc <- paleoTS::as.paleoTSfit(logL = NA, parameters = NA, modelName = paste("OU",
                                                                      "URW", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                        se = NULL)
    return(wc)
  }
  if (hess)
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- paleoTS::as.paleoTSfit(logL = w$value, parameters = w$par, modelName = paste("OU",
                                                                            "URW", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                      se = w$se)
  return(wc)
}


logL.joint.OU.URW<-function (p, y, gg)
{
  #These parameters need to mach the order of p0
  anc<-p[1]
  vstep.OU <- p[2]
  theta.OU <- p[3]
  alpha <- p[4]
  vstep.urw <- p[5]

  n <- length(y$mm)

  #define a vector based on length of first and second RW in the time series
  rw.seg <- which(gg == 2)
  OU.seg <- which(gg == 1)

  #Extract the time vector for the random walks:
  tt.rw <- y$tt[rw.seg] - y$tt[rw.seg[1]]
  #tt.rw <- y$tt[rw.seg]
  tt.OU <- y$tt[OU.seg]

  anc.2<-tail(ou.M(anc, theta.OU, alpha, tt.OU),1)
  #creates a vector that repeates the ancestral value (first sample mean in sequence)
  M <- c(ou.M(anc, theta.OU, alpha, tt.OU), rep(anc.2, length(tt.rw)))
  M <- unname(M)

  #Multiply variance parameter (vstep) from RW with outer product of time*time
  VVrw <- vstep.urw * outer(tt.rw, tt.rw, FUN = pmin)
  ff <- function(a, b) abs(a - b)
  VV.OU <- outer(tt.OU, tt.OU, FUN = ff)
  VV.OU <- exp(-alpha * VV.OU)
  VVd.OU <- ou.V(vstep.OU, alpha, tt.OU)
  VV2.OU <- outer(VVd.OU, VVd.OU, pmin)
  VV.OU <- VV.OU * VV2.OU

  #Create empty n*n matrix
  VVtot <- array(0, dim = c(n, n))

  #Create final variance matrix by combining the variance matrices from stasis and RW
  VVtot[rw.seg, rw.seg] <- VVrw
  VVtot[OU.seg, OU.seg] <- VV.OU

  #Add population variance to the diagonal of the variance matrix
  diag(VVtot) <- diag(VVtot) + y$vv/y$nn

  S <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VVtot, log = TRUE)
  return(S)
}

######

opt.joint.OU.GRW<-function (y, gg, cl = list(fnscale = -1), pool = TRUE, meth = "L-BFGS-B", hess = FALSE)
{
  if (pool) {
    y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  }
  small <- 1e-08
  p0anc<-y$mm[1]
  names(p0anc) <- "anc"
  p0dt <- paleoTS::mle.GRW(paleoTS::sub.paleoTS(y, ok = gg == 2))
  names(p0dt)<-c("mstep.grw", "vstep.grw")
  halft <- (paleoTS::sub.paleoTS(y, ok = gg == 1)$tt[length(paleoTS::sub.paleoTS(y, ok = gg == 1)$tt)] - paleoTS::sub.paleoTS(y, ok = gg == 1)$tt[1])/4
  p0OU <- c(paleoTS::mle.GRW(paleoTS::sub.paleoTS(y, ok = gg == 1))[2]/10, paleoTS::sub.paleoTS(y, ok = gg == 1)$mm[length(paleoTS::sub.paleoTS(y, ok = gg == 1)$mm)], log(2)/halft)
  names(p0OU) <- c("vstep.OU", "theta_OU", "alpha")

  K <- 7

  #Checking if te variance parameter in the RW and Stasis is below critical value, and if TRUE, multiply it with 100
  if (p0dt["vstep.grw"] <= small) p0dt["vstep.grw"] <- 100 * small
  if (p0OU["vstep.OU"] <= small) p0OU["vstep.OU"] <- 100 * small

  #make vector out of estimated model parameters
  p0 <- c(p0anc, p0OU, p0dt)

  #define lower bounds for the L-BFGS-B method
  # three parameters in the unbiased ramdom walk and 4 parameters for the biased random walk
  ll<- c(NA, small, NA, small,NA, small)

  # If user wants to use a differen method than L-BFGS-B, then the lower boun is defined as -Inf
  if (meth != "L-BFGS-B")
    ll <- -Inf

  #Here the multivariate parameter estimation routine start:
  w <- try(optim(p0, fn = logL.joint.OU.GRW, gg = gg,
                 method = meth, lower = ll, control = cl, hessian = hess,
                 y = y), silent = TRUE)
  if (inherits (w, "try-error")) {
    cl <- list(fnscale = -1, parscale = c(1, 10, 100))
    w <- try(optim(p0, fn = logL.joint.OU.GRW, gg = gg,
                   method = meth, lower = ll, control = cl, hessian = hess,
                   y = y), silent = TRUE)
  }

  if (inherits (w, "try-error")) {
    wc <- paleoTS::as.paleoTSfit(logL = NA, parameters = NA, modelName = paste("OU",
                                                                      "GRW", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                        se = NULL)
    return(wc)
  }
  if (hess)
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- paleoTS::as.paleoTSfit(logL = w$value, parameters = w$par, modelName = paste("OU",
                                                                            "GRW", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                      se = w$se)
  return(wc)
}

logL.joint.OU.GRW<-function (p, y, gg)
{
  #These parameters need to mach the order of p0
  anc<-p[1]
  vstep.OU <- p[2]
  theta.OU <- p[3]
  alpha <- p[4]
  mstep.grw <- p[5]
  vstep.grw <- p[6]

  n <- length(y$mm)

  #define a vector based on length of first and second RW in the time series
  dt.seg <- which(gg == 2)
  OU.seg <- which(gg == 1)

  #Extract the time vector for the random walks:
  tt.dt <- y$tt[dt.seg] - y$tt[dt.seg[1]]
  #tt.dt <- y$tt[dt.seg]
  tt.OU <- y$tt[OU.seg]

  anc.2<-tail(ou.M(anc, theta.OU, alpha, tt.OU),1)
  #creates a vector that repeates the ancestral value (first sample mean in sequence)
  M <- c(ou.M(anc, theta.OU, alpha, tt.OU), (anc.2 + mstep.grw * tt.dt))
  M <- unname(M)

  #Multiply variance parameter (vstep) from RW with outer product of time*time
  VV.dt <- vstep.grw * outer(tt.dt, tt.dt, FUN = pmin)
  ff <- function(a, b) abs(a - b)
  VV.OU <- outer(tt.OU, tt.OU, FUN = ff)
  VV.OU <- exp(-alpha * VV.OU)
  VVd.OU <- ou.V(vstep.OU, alpha, tt.OU)
  VV2.OU <- outer(VVd.OU, VVd.OU, pmin)
  VV.OU <- VV.OU * VV2.OU

  #Create empty n*n matrix
  VVtot <- array(0, dim = c(n, n))

  #Create final variance matrix by combining the variance matrices from stasis and RW
  VVtot[dt.seg, dt.seg] <- VV.dt
  VVtot[OU.seg, OU.seg] <- VV.OU

  #Add population variance to the diagonal of the variance matrix
  diag(VVtot) <- diag(VVtot) + y$vv/y$nn

  S <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VVtot, log = TRUE)
  return(S)
}

######


opt.joint.OU.OU<-function (y, gg, cl = list(fnscale = -1), pool = TRUE, meth = "L-BFGS-B", hess = FALSE)
{
  if (pool) {
    y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  }
  small <- 1e-08
  p0OU_1<- paleoTS::mle.GRW(paleoTS::sub.paleoTS(y, ok = gg == 1))
  p0OU_2<- paleoTS::mle.GRW(paleoTS::sub.paleoTS(y, ok = gg == 2))

  K <- 8

  p0anc<-y$mm[1]
  names(p0anc)<-"anc"

  #prepare initial guesses of OU parameters
  halft_1 <- (paleoTS::sub.paleoTS(y, ok = gg == 1)$tt[length(paleoTS::sub.paleoTS(y, ok = gg == 1)$tt)] - paleoTS::sub.paleoTS(y, ok = gg == 1)$tt[1])/4
  p0_1 <- c(p0OU_1["vstep"]/10, paleoTS::sub.paleoTS(y, ok = gg == 1)$mm[length(paleoTS::sub.paleoTS(y, ok = gg == 1)$mm)], log(2)/halft_1)
  names(p0_1) <- c("vstep_1", "theta_1", "alpha_1")

  halft_2 <- (paleoTS::sub.paleoTS(y, ok = gg == 2)$tt[length(paleoTS::sub.paleoTS(y, ok = gg == 2)$tt)] - paleoTS::sub.paleoTS(y, ok = gg == 2)$tt[1])/4
  p0_2 <- c(p0OU_2["vstep"]/10, paleoTS::sub.paleoTS(y, ok = gg == 2)$mm[length(paleoTS::sub.paleoTS(y, ok = gg == 2)$mm)], log(2)/halft_2)
  names(p0_2) <- c("vstep_2", "theta_2", "alpha_2")

  #make vector out of estimated model parameters
  p0 <- c(p0anc, p0_1, p0_2 )

  #define lower bounds for the L-BFGS-B method
  # three parameters in the unbiased ramdom walk and 4 parameters for the biased random walk
  ll<- c(NA, small, NA, small, small, NA, small)

  # If user wants to use a differen method than L-BFGS-B, then the lower boun is defined as -Inf
  if (meth != "L-BFGS-B")
    ll <- -Inf

  #Here the multivariate parameter estimation routine start:
  w <- try(optim(p0, fn = logL.joint.OU.OU, gg = gg,
                 method = meth, lower = ll, control = cl, hessian = hess,
                 y = y), silent = TRUE)
  if (inherits (w, "try-error")) {
    cl <- list(fnscale = -1, parscale = c(1, 10, 100))
    w <- try(optim(p0, fn = logL.joint.OU.OU, gg = gg,
                   method = meth, lower = ll, control = cl, hessian = hess,
                   y = y), silent = TRUE)
  }

  if (inherits (w, "try-error")) {
    wc <- paleoTS::as.paleoTSfit(logL = NA, parameters = NA, modelName = paste("OU",
                                                                      "OU", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                        se = NULL)
    return(wc)
  }
  if (hess)
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- paleoTS::as.paleoTSfit(logL = w$value, parameters = w$par, modelName = paste("OU",
                                                                            "OU", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                      se = w$se)
  return(wc)
}

logL.joint.OU.OU<-function (p, y, gg)
{
  #These parameters need to match the order of p0
  anc<-      p[1]
  vstep_1 <- p[2]
  theta_1 <- p[3]
  alpha_1 <- p[4]
  vstep_2 <- p[5]
  theta_2 <- p[6]
  alpha_2 <- p[7]

  n <- length(y$mm)

  #define a vector based on length of first and second RW in the time series
  OU.seg_1 <- which(gg == 1)
  OU.seg_2 <- which(gg == 2)

  #Extract the time vector for the OU:
  tt.OU_1 <- y$tt[OU.seg_1]
  tt.OU_2 <- y$tt[OU.seg_2] - y$tt[OU.seg_2[1]]

  #create a vector for the expected means for the whole time series
  M <- c(ou.M(anc, theta_1, alpha_1, tt.OU_1), ou.M(theta_1, theta_2, alpha_2, tt.OU_2))
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

  #Create empty n*n matrix
  VVtot <- array(0, dim = c(n, n))

  #Create final variance matrix by combining the variance matrices from stasis and RW
  VVtot[OU.seg_1, OU.seg_1] <- VV.OU_1
  VVtot[OU.seg_2, OU.seg_2] <- VV.OU_2

  #Add population variance to the diagonal of the variance matrix
  diag(VVtot) <- diag(VVtot) + y$vv/y$nn

  #
  S <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VVtot, log = TRUE)
  return(S)
}

###########

opt.joint.GRW.OU<-function (y, gg, cl = list(fnscale = -1), pool = TRUE, meth = "L-BFGS-B", hess = FALSE)
{
  if (pool) {
    y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  }
  small <- 1e-08
  p0anc<-y$mm[1]
  names(p0anc) <- "anc"
  p0dt <- paleoTS::mle.GRW(paleoTS::sub.paleoTS(y, ok = gg == 1))
  names(p0dt)<-c("mstep.grw", "vstep.grw")
  halft <- (paleoTS::sub.paleoTS(y, ok = gg == 2)$tt[length(paleoTS::sub.paleoTS(y, ok = gg == 2)$tt)] - paleoTS::sub.paleoTS(y, ok = gg == 2)$tt[1])/4
  p0OU <- c(paleoTS::mle.GRW(paleoTS::sub.paleoTS(y, ok = gg == 2))[2]/10, paleoTS::sub.paleoTS(y, ok = gg == 2)$mm[length(paleoTS::sub.paleoTS(y, ok = gg == 2)$mm)], log(2)/halft)
  names(p0OU) <- c("vstep.OU", "theta_OU", "alpha")

  K <- 7

  #Checking if te variance parameter in the RW and Stasis is below critical value, and if TRUE, multiply it with 100
  if (p0dt["vstep.grw"] <= small) p0dt["vstep.grw"] <- 100 * small
  if (p0OU["vstep.OU"] <= small) p0OU["vstep.OU"] <- 100 * small

  #make vector out of estimated model parameters
  p0 <- c(p0anc, p0dt, p0OU)

  #define lower bounds for the L-BFGS-B method
  # three parameters in the unbiased ramdom walk and 4 parameters for the biased random walk
  ll<- c(NA, NA, small, small,NA, small)

  # If user wants to use a differen method than L-BFGS-B, then the lower boun is defined as -Inf
  if (meth != "L-BFGS-B")
    ll <- -Inf

  #Here the multivariate parameter estimation routine start:
  w <- try(optim(p0, fn = logL.joint.GRW.OU, gg = gg,
                 method = meth, lower = ll, control = cl, hessian = hess,
                 y = y), silent = TRUE)
  if (inherits (w, "try-error")) {
    cl <- list(fnscale = -1, parscale = c(1, 10, 100))
    w <- try(optim(p0, fn = logL.joint.GRW.OU, gg = gg,
                   method = meth, lower = ll, control = cl, hessian = hess,
                   y = y), silent = TRUE)
  }

  if (inherits (w, "try-error")) {
    wc <- paleoTS::as.paleoTSfit(logL = NA, parameters = NA, modelName = paste("GRW",
                                                                      "OU", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                        se = NULL)
    return(wc)
  }
  if (hess)
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- paleoTS::as.paleoTSfit(logL = w$value, parameters = w$par, modelName = paste("GRW",
                                                                            "OU", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                      se = w$se)
  return(wc)
}

logL.joint.GRW.OU<-function (p, y, gg)
{
  #These parameters need to mach the order of p0
  anc<-p[1]
  mstep.grw <- p[2]
  vstep.grw <- p[3]
  vstep.OU <- p[4]
  theta.OU <- p[5]
  alpha <- p[6]

  n <- length(y$mm)

  #define a vector based on length of first and second RW in the time series
  dt.seg <- which(gg == 1)
  OU.seg <- which(gg == 2)

  #Extract the time vector for the random walks:
  tt.dt <- y$tt[dt.seg]
  tt.OU <- y$tt[OU.seg] - y$tt[OU.seg[1]]
  #tt.OU <- y$tt[OU.seg]

  anc.2<-tail((anc + mstep.grw * tt.dt),1)
  #creates a vector that repeates the ancestral value (first sample mean in sequence)
  M <- c((anc + mstep.grw * tt.dt), ou.M(anc.2, theta.OU, alpha, tt.OU))
  M <- unname(M)

  #Multiply variance parameter (vstep) from RW with outer product of time*time
  VVdt <- vstep.grw * outer(tt.dt, tt.dt, FUN = pmin)
  ff <- function(a, b) abs(a - b)
  VV.OU <- outer(tt.OU, tt.OU, FUN = ff)
  VV.OU <- exp(-alpha * VV.OU)
  VVd.OU <- ou.V(vstep.OU, alpha, tt.OU)
  VV2.OU <- outer(VVd.OU, VVd.OU, pmin)
  VV.OU <- VV.OU * VV2.OU

  #Create empty n*n matrix
  VVtot <- array(0, dim = c(n, n))

  #Create final variance matrix by combining the variance matrices from stasis and RW
  VVtot[dt.seg, dt.seg] <- VVdt
  VVtot[OU.seg, OU.seg] <- VV.OU

  #Add population variance to the diagonal of the variance matrix
  diag(VVtot) <- diag(VVtot) + y$vv/y$nn

  S <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VVtot, log = TRUE)
  return(S)
}

#############

opt.joint.GRW.GRW<-function (y, gg, cl = list(fnscale = -1), pool = TRUE, meth = "L-BFGS-B", hess = FALSE)
{
  if (pool) {
    y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  }
  small <- 1e-08
  p0dt.1 <- paleoTS::mle.GRW(paleoTS::sub.paleoTS(y, ok = gg == 1))
  p0dt.2 <- paleoTS::mle.GRW(paleoTS::sub.paleoTS(y, ok = gg == 2))

  K <- 6

  #Checking if te variance parameter in the RW and Stasis is below critical value, and if TRUE, multiply it with 100
  if (p0dt.1["vstep"] <= small)
    p0dt.1["vstep"] <- 100 * small
  if (p0dt.2["vstep"] <= small)
    p0dt.2["vstep"] <- 100 * small

  # Define ancestral state
  p0anc <- y$mm[1]
  names(p0anc)  <- "anc"
  names(p0dt.1) <- c("mstep.1", "vstep.1")
  names(p0dt.2) <- c("mstep.2", "vstep.2")

  #make vector out of estimated model parameters
  p0 <- c(p0anc, p0dt.1, p0dt.2)

  #define lower bounds for the L-BFGS-B method
  # three parameters in the unbiased ramdom walk and 4 parameters for the biased random walk
  ll<- c(NA, NA, small, NA, small)

  # If user wants to use a differen method than L-BFGS-B, then the lower boun is defined as -Inf
  if (meth != "L-BFGS-B")
    ll <- -Inf

  #Here the multivariate parameter estimation routine start:
  w <- try(optim(p0, fn = logL.joint.GRW.GRW, gg = gg,
                 method = meth, lower = ll, control = cl, hessian = hess,
                 y = y), silent = TRUE)
  if (inherits (w, "try-error")) {
    cl <- list(fnscale = -1, parscale = c(1, 10, 100))
    w <- try(optim(p0, fn = logL.joint.GRW.GRW, gg = gg,
                   method = meth, lower = ll, control = cl, hessian = hess,
                   y = y), silent = TRUE)
  }

  if (inherits (w, "try-error")) {
    wc <- paleoTS::as.paleoTSfit(logL = NA, parameters = NA, modelName = paste("GRW",
                                                                      "GRW", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                        se = NULL)
    return(wc)
  }
  if (hess)
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- paleoTS::as.paleoTSfit(logL = w$value, parameters = w$par, modelName = paste("GRW",
                                                                            "GRW", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                      se = w$se)
  return(wc)
}

logL.joint.GRW.GRW<-function (p, y, gg)
{
  #These parameters need to mach the order of p0
  anc<-p[1]
  mstep.1 <- p[2]
  vstep.1 <- p[3]
  mstep.2 <- p[4]
  vstep.2 <- p[5]

  n <- length(y$mm)

  #define a vector based on length of first and second RW in the time series
  dt.1.seg <- which(gg == 1)
  dt.2.seg <- which(gg == 2)

  #Extract the time vector for the random walks:
  tt.dt.1 <- y$tt[dt.1.seg]
  tt.dt.2 <- y$tt[dt.2.seg] - y$tt[dt.2.seg[1]]

  anc.2<-tail((anc + mstep.1 * tt.dt.1),1)

  #creates a vector that repeates the ancestral value (first sample mean in sequence)
  M <- c((anc + mstep.1 * tt.dt.1), (anc.2 + mstep.2 * tt.dt.2))
  M <- unname(M)

  #Multiply variance parameter (vstep) from RW with outer product of time*time
  VVdt.1 <- vstep.1 * outer(tt.dt.1, tt.dt.1, FUN = pmin)
  VVdt.2 <- vstep.2 * outer(tt.dt.2, tt.dt.2, FUN = pmin)

  #Create empty n*n matrix
  VVtot <- array(0, dim = c(n, n))

  #Create final variance matrix by combining the variance matrices from stasis and RW
  VVtot[dt.1.seg, dt.1.seg] <- VVdt.1
  VVtot[dt.2.seg, dt.2.seg] <- VVdt.2

  #Add population variance to the diagonal of the variance matrix
  diag(VVtot) <- diag(VVtot) + y$vv/y$nn

  S <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VVtot, log = TRUE)
  return(S)
}

###########

opt.joint.URW.OU<-function (y, gg, cl = list(fnscale = -1), pool = TRUE, meth = "L-BFGS-B", hess = FALSE)
{
  if (pool) {
    y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  }
  small <- 1e-08
  p0anc <- y$mm[1]
  p0rw <- paleoTS::mle.URW(paleoTS::sub.paleoTS(y, ok = gg == 1))
  names(p0anc) <- "anc"
  names(p0rw)<-"vstep.urw"
  halft <- (paleoTS::sub.paleoTS(y, ok = gg == 2)$tt[length(paleoTS::sub.paleoTS(y, ok = gg == 2)$tt)] - paleoTS::sub.paleoTS(y, ok = gg == 2)$tt[1])/4
  p0OU <- c(paleoTS::mle.GRW(paleoTS::sub.paleoTS(y, ok = gg == 2))[2]/10, paleoTS::sub.paleoTS(y, ok = gg == 2)$mm[length(paleoTS::sub.paleoTS(y, ok = gg == 2)$mm)], log(2)/halft)
  names(p0OU) <- c("vstep.OU", "theta_OU", "alpha")

  K <- 6

  #Checking if te variance parameter in the RW and Stasis is below critical value, and if TRUE, multiply it with 100
  if (p0rw["vstep.urw"] <= small) p0rw["vstep.urw"] <- 100 * small
  if (p0OU["vstep.OU"] <= small) p0OU["vstep.OU"] <- 100 * small

  #make vector out of estimated model parameters
  p0 <- c(p0anc, p0rw, p0OU)

  #define lower bounds for the L-BFGS-B method
  # three parameters in the unbiased ramdom walk and 4 parameters for the biased random walk
  ll<- c(NA, small, small,NA, small)

  # If user wants to use a differen method than L-BFGS-B, then the lower boun is defined as -Inf
  if (meth != "L-BFGS-B")
    ll <- -Inf

  #Here the multivariate parameter estimation routine start:
  w <- try(optim(p0, fn = logL.joint.URW.OU, gg = gg,
                 method = meth, lower = ll, control = cl, hessian = hess,
                 y = y), silent = TRUE)
  if (inherits (w, "try-error")) {
    cl <- list(fnscale = -1, parscale = c(1, 10, 100))
    w <- try(optim(p0, fn = logL.joint.URW.OU, gg = gg,
                   method = meth, lower = ll, control = cl, hessian = hess,
                   y = y), silent = TRUE)
  }

  if (inherits (w, "try-error")) {
    wc <- paleoTS::as.paleoTSfit(logL = NA, parameters = NA, modelName = paste("URW",
                                                                      "OU", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                        se = NULL)
    return(wc)
  }
  if (hess)
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- paleoTS::as.paleoTSfit(logL = w$value, parameters = w$par, modelName = paste("URW",
                                                                            "OU", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                      se = w$se)
  return(wc)
}

logL.joint.URW.OU<-function (p, y, gg)
{
  #These parameters need to mach the order of p0
  anc<-p[1]
  vstep.urw <- p[2]
  vstep.OU <- p[3]
  theta.OU <- p[4]
  alpha <- p[5]

  n <- length(y$mm)

  #define a vector based on length of first and second RW in the time series
  rw.seg <- which(gg == 1)
  OU.seg <- which(gg == 2)

  #Extract the time vector for the random walks:
  tt.rw <- y$tt[rw.seg]
  tt.OU <- y$tt[OU.seg] - y$tt[OU.seg[1]]
  #tt.OU <- y$tt[OU.seg]

  #creates a vector that repeates the ancestral value (first sample mean in sequence)
  M <- c(rep(anc, length(tt.rw)), ou.M(anc, theta.OU, alpha, tt.OU))
  M <- unname(M)

  #Multiply variance parameter (vstep) from RW with outer product of time*time
  VVrw <- vstep.urw * outer(tt.rw, tt.rw, FUN = pmin)
  ff <- function(a, b) abs(a - b)
  VV.OU <- outer(tt.OU, tt.OU, FUN = ff)
  VV.OU <- exp(-alpha * VV.OU)
  VVd.OU <- ou.V(vstep.OU, alpha, tt.OU)
  VV2.OU <- outer(VVd.OU, VVd.OU, pmin)
  VV.OU <- VV.OU * VV2.OU

  #Create empty n*n matrix
  VVtot <- array(0, dim = c(n, n))

  #Create final variance matrix by combining the variance matrices from stasis and RW
  VVtot[rw.seg, rw.seg] <- VVrw
  VVtot[OU.seg, OU.seg] <- VV.OU

  #Add population variance to the diagonal of the variance matrix
  diag(VVtot) <- diag(VVtot) + y$vv/y$nn

  S <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VVtot, log = TRUE)
  return(S)
}

#############

opt.joint.URW.GRW<-function (y, gg, cl = list(fnscale = -1), pool = TRUE, meth = "L-BFGS-B", hess = FALSE)
{
  if (pool) {
    y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  }
  small <- 1e-08
  p0rw <- paleoTS::mle.URW(paleoTS::sub.paleoTS(y, ok = gg == 1))
  p0dt <- paleoTS::mle.GRW(paleoTS::sub.paleoTS(y, ok = gg == 2))

  K <- 5

  #Checking if te variance parameter in the RW and Stasis is below critical value, and if TRUE, multiply it with 100
  if (p0rw["vstep"] <= small)
    p0rw["vstep"] <- 100 * small
  if (p0dt["vstep"] <= small)
    p0dt["vstep"] <- 100 * small

  # Define ancestral state
  p0anc <- y$mm[1]
  names(p0anc) <- "anc"
  names(p0rw)<-"vstep.urw"
  names(p0dt)<-c("mstep", "vstep.grw")

  #make vector out of estimated model parameters
  p0 <- c(p0anc, p0rw, p0dt)

  #define lower bounds for the L-BFGS-B method
  # three parameters in the unbiased ramdom walk and 4 parameters for the biased random walk
  ll<- c(NA, small, NA, small)

  # If user wants to use a differen method than L-BFGS-B, then the lower boun is defined as -Inf
  if (meth != "L-BFGS-B")
    ll <- -Inf

  #Here the multivariate parameter estimation routine start:
  w <- try(optim(p0, fn = logL.joint.URW.GRW, gg = gg,
                 method = meth, lower = ll, control = cl, hessian = hess,
                 y = y), silent = TRUE)
  if (inherits (w, "try-error")) {
    cl <- list(fnscale = -1, parscale = c(1, 10, 100))
    w <- try(optim(p0, fn = logL.joint.URW.GRW, gg = gg,
                   method = meth, lower = ll, control = cl, hessian = hess,
                   y = y), silent = TRUE)
  }

  if (inherits (w, "try-error")) {
    wc <- paleoTS::as.paleoTSfit(logL = NA, parameters = NA, modelName = paste("URW",
                                                                      "GRW", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                        se = NULL)
    return(wc)
  }
  if (hess)
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- paleoTS::as.paleoTSfit(logL = w$value, parameters = w$par, modelName = paste("URW",
                                                                            "GRW", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                      se = w$se)
  return(wc)
}

logL.joint.URW.GRW<-function (p, y, gg)
{
  #These parameters need to mach the order of p0
  anc<-p[1]
  vstep.urw <- p[2]
  mstep <- p[3]
  vstep.grw <- p[4]

  n <- length(y$mm)

  #define a vector based on length of first and second RW in the time series
  rw.seg <- which(gg == 1)
  dt.seg <- which(gg == 2)

  #Extract the time vector for the random walks:
  tt.rw <- y$tt[rw.seg]
  tt.dt <- y$tt[dt.seg] - y$tt[dt.seg[1]]
 #tt.dt <- y$tt[dt.seg]

  #creates a vector that repeates the ancestral value (first sample mean in sequence)
  M <- c(rep(anc, length(tt.rw)), anc + mstep * tt.dt)
  M <- unname(M)

  #Multiply variance parameter (vstep) from RW with outer product of time*time
  VVrw <- vstep.urw * outer(tt.rw, tt.rw, FUN = pmin)
  VVdt <- vstep.grw * outer(tt.dt, tt.dt, FUN = pmin)

  #Create empty n*n matrix
  VVtot <- array(0, dim = c(n, n))

  #Create final variance matrix by combining the variance matrices from stasis and RW
  VVtot[rw.seg, rw.seg] <- VVrw
  VVtot[dt.seg, dt.seg] <- VVdt

  #Add population variance to the diagonal of the variance matrix
  diag(VVtot) <- diag(VVtot) + y$vv/y$nn

  #
  S <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VVtot, log = TRUE)
  return(S)
}

#############

opt.joint.Stasis.GRW<-function (y, gg, cl = list(fnscale = -1), pool = TRUE, meth = "L-BFGS-B", hess = FALSE)
{
  if (pool) {
    y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  }
  small <- 1e-08
  p0st <- paleoTS::mle.Stasis(paleoTS::sub.paleoTS(y, ok = gg == 1))
  p0dt<- paleoTS::mle.GRW(paleoTS::sub.paleoTS(y, ok = gg == 2))

  K <- 5

  #Checking if te variance parameter in the RW and Stasis is below critical value, and if TRUE, multiply it with 100
  if (p0st["omega"] <= small) p0st["omega"] <- 100 * small
  if (p0dt["vstep"] <= small) p0dt["vstep"] <- 100 * small

  #make vector out of estimated model parameters
  p0 <- c(p0st, p0dt)

  #define lower bounds for the L-BFGS-B method
  # three parameters in the unbiased ramdom walk and 4 parameters for the biased random walk
  ll<- c(NA, small, NA, small)

  # If user wants to use a differen method than L-BFGS-B, then the lower bound is defined as -Inf
  if (meth != "L-BFGS-B")
    ll <- -Inf

  #Here the multivariate parameter estimation routine start:
  w <- try(optim(p0, fn = logL.joint.Stasis.GRW, gg = gg,
                 method = meth, lower = ll, control = cl, hessian = hess,
                 y = y), silent = TRUE)
  if (inherits (w, "try-error")) {
    cl <- list(fnscale = -1, parscale = c(1, 10, 100))
    w <- try(optim(p0, fn = logL.joint.Stasis.GRW, gg = gg,
                   method = meth, lower = ll, control = cl, hessian = hess,
                   y = y), silent = TRUE)
  }

  if (inherits (w, "try-error")) {
    wc <- paleoTS::as.paleoTSfit(logL = NA, parameters = NA, modelName = paste("Stasis",
                                                                      "GRW", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                        se = NULL)
    return(wc)
  }
  if (hess)
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- paleoTS::as.paleoTSfit(logL = w$value, parameters = w$par, modelName = paste("Stasis",
                                                                            "GRW", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                      se = w$se)
  return(wc)
}

logL.joint.Stasis.GRW<-function (p, y, gg)
{
  #These parameters need to match the order of p0
  theta <- p[1]
  omega <- p[2]
  mstep <- p[3]
  vstep <- p[4]

  n <- length(y$mm)

  #define a vector based on length of first and second RW in the time series
  st.seg <- which(gg == 1)
  dt.seg <- which(gg == 2)

  #Extract the time vector for the trend:
  #tt.trend <- y$tt[dt.seg]
  tt.trend <- y$tt[dt.seg] - y$tt[dt.seg[1]]

  #create a vector for the expected means for the whole time series
  M <- c(rep(theta, length(st.seg)), theta + mstep * tt.trend)
  M <- unname(M)

  #Create variance-covariance matrices for the three models
  VVst <- diag(omega, nrow = length(st.seg))
  VVtrend <- vstep * outer(tt.trend, tt.trend, FUN = pmin)

  #Create empty n*n matrix
  VVtot <- array(0, dim = c(n, n))

  #Create final variance matrix by combining the variance matrices from stasis and RW
  VVtot[st.seg, st.seg] <- VVst
  VVtot[dt.seg, dt.seg] <- VVtrend

  #Add population variance to the diagonal of the variance matrix
  diag(VVtot) <- diag(VVtot) + y$vv/y$nn

  S <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VVtot, log = TRUE)
  return(S)
}


############

opt.joint.Stasis.URW<-function (y, gg, cl = list(fnscale = -1), pool = TRUE, meth = "L-BFGS-B", hess = FALSE)
{
  if (pool) {
    y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  }
  small <- 1e-08
  p0st <- paleoTS::mle.Stasis(paleoTS::sub.paleoTS(y, ok = gg == 1))
  names(p0st) <- c("theta", "omega")
  p0urw <- paleoTS::mle.URW(paleoTS::sub.paleoTS(y, ok = gg == 2))
  names(p0urw) <-"vstep"

  K <- 4

  #Checking if te variance parameter in the RW and Stasis is below critical value, and if TRUE, multiply it with 100
  if (p0st["omega"] <= small) p0st["omega"] <- 100 * small
  if (p0urw["vstep"] <= small) p0urw["vstep"] <- 100 * small

  #make vector out of estimated model parameters
  p0 <- c(p0st, p0urw)

  #define lower bounds for the L-BFGS-B method
  # three parameters in the unbiased ramdom walk and 4 parameters for the biased random walk
  ll<- c(NA, small, small)

  # If user wants to use a differen method than L-BFGS-B, then the lower boun is defined as -Inf
  if (meth != "L-BFGS-B")
    ll <- -Inf

  #Here the multivariate parameter estimation routine start:
  w <- try(optim(p0, fn = logL.joint.Stasis.URW, gg = gg,
                 method = meth, lower = ll, control = cl, hessian = hess,
                 y = y), silent = TRUE)
  if (inherits (w, "try-error")) {
    cl <- list(fnscale = -1, parscale = c(1, 10, 100))
    w <- try(optim(p0, fn = logL.joint.Stasis.URW, gg = gg,
                   method = meth, lower = ll, control = cl, hessian = hess,
                   y = y), silent = TRUE)
  }

  if (inherits (w, "try-error")) {
    wc <- paleoTS::as.paleoTSfit(logL = NA, parameters = NA, modelName = paste("Stasis",
                                                                      "URW", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                        se = NULL)
    return(wc)
  }
  if (hess)
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- paleoTS::as.paleoTSfit(logL = w$value, parameters = w$par, modelName = paste("Stasis",
                                                                            "URW", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                      se = w$se)
  return(wc)
}

logL.joint.Stasis.URW<-function (p, y, gg)
{
  theta <- p[1]
  omega <- p[2]
  vstep <- p[3]

  n <- length(y$mm)
  M <- rep(theta, n)

  #define a vector based on length of first and second RW in the time series
  st.seg <- which(gg == 1)
  urw.seg <- which(gg == 2)
  tt.urw <- y$tt[urw.seg] - y$tt[urw.seg[1]]
  #tt.urw <- y$tt[urw.seg]

  #Create variance-covariance matrices for the three models
  VV.st <- diag(omega, nrow = length(st.seg))
  VV.urw <- vstep * outer(tt.urw, tt.urw, FUN = pmin)

  #Create empty n*n matrix
  VVtot <- array(0, dim = c(n, n))

  #Create final variance matrix by combining the variance matrices from stasis and RW
  VVtot[st.seg, st.seg] <- VV.st
  VVtot[urw.seg, urw.seg] <- VV.urw

  #Add population variance to the diagonal of the variance matrix
  diag(VVtot) <- diag(VVtot) + y$vv/y$nn

  #
  S <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VVtot, log = TRUE)
  return(S)
}


#######
opt.joint.Stasis.Stasis<-function (y, gg, cl = list(fnscale = -1), pool = TRUE, meth = "L-BFGS-B", hess = FALSE)
{
  if (pool) {
    y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  }
  small <- 1e-08
  p0st.1 <- paleoTS::mle.Stasis(paleoTS::sub.paleoTS(y, ok = gg == 1))
  names(p0st.1) <- c("theta", "omega.1")
  p0st.2 <- paleoTS::mle.Stasis(paleoTS::sub.paleoTS(y, ok = gg == 2))[2]
  names(p0st.2) <-"omega.2"

  K <- 4

  #Checking if te variance parameter in the RW and Stasis is below critical value, and if TRUE, multiply it with 100
  if (p0st.1["omega.1"] <= small) p0st.1["omega.1"] <- 100 * small
  if (p0st.2["omega.2"] <= small) p0st.2["omega.2"] <- 100 * small

  #make vector out of estimated model parameters
  p0 <- c(p0st.1, p0st.2)

  #define lower bounds for the L-BFGS-B method
  # three parameters in the unbiased ramdom walk and 4 parameters for the biased random walk
  ll<- c(NA, small, small)

  # If user wants to use a differen method than L-BFGS-B, then the lower boun is defined as -Inf
  if (meth != "L-BFGS-B")
    ll <- -Inf

  #Here the multivariate parameter estimation routine start:
  w <- try(optim(p0, fn = logL.joint.Stasis.Stasis, gg = gg,
                 method = meth, lower = ll, control = cl, hessian = hess,
                 y = y), silent = TRUE)
  if (inherits (w, "try-error")) {
    cl <- list(fnscale = -1, parscale = c(1, 10, 100))
    w <- try(optim(p0, fn = logL.joint.Stasis.Stasis, gg = gg,
                   method = meth, lower = ll, control = cl, hessian = hess,
                   y = y), silent = TRUE)
  }

  if (inherits (w, "try-error")) {
    wc <- paleoTS::as.paleoTSfit(logL = NA, parameters = NA, modelName = paste("Stasis",
                                                                      "Stasis", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                        se = NULL)
    return(wc)
  }
  if (hess)
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- paleoTS::as.paleoTSfit(logL = w$value, parameters = w$par, modelName = paste("Stasis",
                                                                            "Stasis", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                      se = w$se)
  return(wc)
}

logL.joint.Stasis.Stasis<-function (p, y, gg)
{
  #These parameters need to match the order of p0
  theta_st <- p[1]
  omega.1 <- p[2]
  omega.2 <- p[3]

  n <- length(y$mm)

  M <- rep(theta_st, n)
  M <- unname(M)

  #define a vector based on length of first and second RW in the time series
  st.1.seg <- which(gg == 1)
  st.2.seg <- which(gg == 2)

  #Create variance-covariance matrices for the three models
  VV.st.1 <- diag(omega.1, nrow = length(st.1.seg))
  VV.st.2 <- diag(omega.2, nrow = length(st.2.seg))

  #Create empty n*n matrix
  VVtot <- array(0, dim = c(n, n))

  #Create final variance matrix by combining the variance matrices from stasis and RW
  VVtot[st.1.seg, st.1.seg] <- VV.st.1
  VVtot[st.2.seg, st.2.seg] <- VV.st.2

  #Add population variance to the diagonal of the variance matrix
  diag(VVtot) <- diag(VVtot) + y$vv/y$nn

  #
  S <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VVtot, log = TRUE)
  return(S)
}



opt.joint.URW.URW<-function (y, gg, cl = list(fnscale = -1), pool = TRUE, meth = "L-BFGS-B", hess = FALSE)
{
  if (pool) {
    y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  }
  small <- 1e-08
  p0rw_1 <- paleoTS::mle.URW(paleoTS::sub.paleoTS(y, ok = gg == 1))
  p0rw_2 <- paleoTS::mle.URW(paleoTS::sub.paleoTS(y, ok = gg == 2))

  K <- 4

  #Checking if te variance parameter in the RW and Stasis is below critical value, and if TRUE, multiply it with 100
  if (p0rw_1["vstep"] <= small)
    p0rw_1["vstep"] <- 100 * small
  if (p0rw_2["vstep"] <= small)
    p0rw_2["vstep"] <- 100 * small

  # Define ancestral state
  p0anc <- y$mm[1]
  names(p0anc) <- "anc"

  #make vector out of estimated model parameters
  p0 <- c(p0anc, p0rw_1, p0rw_2)

  #define lower bounds for the L-BFGS-B method
  # three parameters in the unbiased ramdom walk and 4 parameters for the biased random walk
  ll<- c(NA, small, small)

  # If user wants to use a differen method than L-BFGS-B, then the lower boun is defined as -Inf
  if (meth != "L-BFGS-B")
    ll <- -Inf

  #Here the multivariate parameter estimation routine start:
  w <- try(optim(p0, fn = logL.joint.URW.URW, gg = gg,
                 method = meth, lower = ll, control = cl, hessian = hess,
                 y = y), silent = TRUE)
  if (inherits (w, "try-error")) {
    cl <- list(fnscale = -1, parscale = c(1, 10, 100))
    w <- try(optim(p0, fn = logL.joint.URW.URW, gg = gg,
                   method = meth, lower = ll, control = cl, hessian = hess,
                   y = y), silent = TRUE)
  }

  if (inherits (w, "try-error")) {
    wc <- paleoTS::as.paleoTSfit(logL = NA, parameters = NA, modelName = paste("URW",
                                                                      "URW", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                        se = NULL)
    return(wc)
  }
  if (hess)
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- paleoTS::as.paleoTSfit(logL = w$value, parameters = w$par, modelName = paste("URW",
                                                                            "URW", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                      se = w$se)
  return(wc)
}


opt.joint.Stasis.OU<-function (y, gg, cl = list(fnscale = -1), pool = TRUE, meth = "L-BFGS-B", hess = FALSE)
{
  if (pool) {
    y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  }
  small <- 1e-08
  p0st <- paleoTS::mle.Stasis(paleoTS::sub.paleoTS(y, ok = gg == 1))
  names(p0st) <- c("theta_st", "omega")
  p0OU<- paleoTS::mle.GRW(paleoTS::sub.paleoTS(y, ok = gg == 2))

  K <- 6

  #Checking if te variance parameter in the RW and Stasis is below critical value, and if TRUE, multiply it with 100
  if (p0st["omega"] <= small)
    p0st["omega"] <- 100 * small
  if (p0OU["vstep"] <= small)
    p0OU["vstep"] <- 100 * small

  #prepare initial guesses of OU parameters
  halft <- (paleoTS::sub.paleoTS(y, ok = gg == 2)$tt[length(paleoTS::sub.paleoTS(y, ok = gg == 2)$tt)] - paleoTS::sub.paleoTS(y, ok = gg == 2)$tt[1])/4
  p0 <- c(p0OU["vstep"]/10, paleoTS::sub.paleoTS(y, ok = gg == 2)$mm[length(paleoTS::sub.paleoTS(y, ok = gg == 2)$mm)], log(2)/halft)
  names(p0) <- c("vstep", "theta_OU", "alpha")

  #make vector out of estimated model parameters
  p0 <- c(p0st, p0)

  #define lower bounds for the L-BFGS-B method
  # three parameters in the unbiased ramdom walk and 4 parameters for the biased random walk
  ll<- c(NA, small, small, NA, small)

  # If user wants to use a differen method than L-BFGS-B, then the lower boun is defined as -Inf
  if (meth != "L-BFGS-B")
    ll <- -Inf

  #Here the multivariate parameter estimation routine start:
  w <- try(optim(p0, fn = logL.joint.Stasis.OU, gg = gg,
                 method = meth, lower = ll, control = cl, hessian = hess,
                 y = y), silent = TRUE)
  if (inherits (w, "try-error")) {
    cl <- list(fnscale = -1, parscale = c(1, 10, 100))
    w <- try(optim(p0, fn = logL.joint.Stasis.OU, gg = gg,
                   method = meth, lower = ll, control = cl, hessian = hess,
                   y = y), silent = TRUE)
  }

  if (inherits (w, "try-error")) {
    wc <- paleoTS::as.paleoTSfit(logL = NA, parameters = NA, modelName = paste("Stasis",
                                                                      "OU", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                        se = NULL)
    return(wc)
  }
  if (hess)
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- paleoTS::as.paleoTSfit(logL = w$value, parameters = w$par, modelName = paste("Stasis",
                                                                            "OU", sep = "-"), method = "Joint", K = K, n = length(y$mm),
                      se = w$se)
  return(wc)
}


#################

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
  tt.OU <- y$tt[OU.seg] - y$tt[OU.seg[1]]

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


###########

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
