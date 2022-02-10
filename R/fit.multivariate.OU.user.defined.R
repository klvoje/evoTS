#' @title Fit user-defined multivariate Ornstein-Uhlenbeck models to multivariate evolutionary sequence (time-series) data.
#'
#' @description Function to find maximum likelihood solutions to a multivariate Ornstein-Uhlenbeck model fitted using user-defined A and R matrices.
#'
#' @param yy a multivariate evoTS object.
#'
#' @param A.user the pull matrix. A user-defined A matrix.
#'
#' @param R.user the drift matrix. A user-defined R matrix.
#'
#' @param method optimization method, passed to function optim. Default is "Nelder-Mead".
#'
#' @param hess logical, indicating whether to calculate standard errors from the Hessian matrix.
#'
#' @param trace logical, indicating whether information on the progress of the optimization is printed.
#'
#' @param iterations the number of times the optimization method is run from different starting points. Default is NULL, meaning the optimization is run once.
#'
#' @param iter.sd defines the standard deviation of the Gaussian distribution from which starting values for the optimization routine is run. Default is 1.
#'
#' @details This function provides users the flexibility to define their own A and R matrices. The possibility to define any A matrices enable detailed investigation of specific evolutionary hypotheses. The parameters to be estimated in the matrices are indicated by the value 1. All other entries in the matrix must be 0.
#'
#' The function searches - using an optimization routine - for the maximum-likelihood solution for the chosen multivariate Ornstein-Uhlenbeck model. The argument 'method' is passed to the 'optim' function and is included for the convenience of users to better control the optimization routine. Note that the the default method (Nelder-Mead) seems to work for most evolutionary sequences. The method L-BFGS-B allows box-constraints on some parameters (e.g. non-negative variance parameters) and is faster than Nelder-Mead, but is less stable than the default method (Nelder-Mead).
#'
#' Initial estimates to start the optimization come from maximum-likelihood estimates of the univariate Ornstein-Uhlenbeck model (from the paleoTS package) fitted to each time-series separately.
#'
#' It is good practice to repeat any numerical optimization procedure from different starting points. This is especially important for complex models as the log-likelihood surface might contain more than one peak. The number of iterations is controlled by the argument 'iterations'. The function will report the model parameters from the iteration with the highest log-likelihood.
#'
#'@return First part of the output reports the log-likelihood of the model and its AICc score. The second part of the output is the maximum log-likelihood model parameters (ancestral.values, optima, A, and R). The half-life is also provided, which is the  The last part of the output gives information about the number of parameters in the model (K), number of samples in the data (n) and number of times the optimization routine was run (iter).
#'
#'@note The models have been implemented to be compatible with the joint parameterization routine in the package paleoTS. The optimization is therefore fit using the actual sample values, with the autocorrelation among samples accounted for in the log-likelihood function. The joint distribution of sample means is multivariate normal, with means and variance-covariances determined by evolutionary parameters and sampling errors.
#'
#'@author Kjetil Lysne Voje
#'
#'@references Reitan, T., Schweder, T. & Henderiks, J. Phenotypic evolution studied by layered stochastic differential equations. \emph{Ann Appl Statistics} 6, 1531–1551 (2012).
#'@references Bartoszek, K., Pienaar, J., Mostad, P., Andersson, S. & Hansen, T. F. A phylogenetic comparative method for studying multivariate adaptation. \emph{J Theor Biol} 314, 204–215 (2012).
#'@references Clavel, J., Escarguel, G. & Merceron, G. mvmorph: an r package for fitting multivariate evolutionary models to morphometric data. \emph{Methods Ecol Evol 6}, 1311–1319 (2015).
#'
#'@export
#'
#'@examples
#'## Generate a evoTS objects by simulating a multivariate dataset
#'x <- sim.multi.OU(15)
#'
#'## Define an A matrix that is lower diagonal.
#'A <- matrix(c(1,0,1,1), nrow=2, byrow=TRUE)
#'
#'## Define a diagonal R matrix.
#'R <- matrix(c(1,0,0,1), nrow=2, byrow=TRUE)
#'
#'## Fit the multivariate Ornstein-Uhlenbeck model to the data.
#'##fit.multivariate.OU.user.defined(x, A.user=A, R.user=R, trace=TRUE)
#'

fit.multivariate.OU.user.defined<-function (yy, A.user=NULL, R.user=NULL, method="Nelder-Mead", hess = FALSE, trace=FALSE, iterations=NULL, iter.sd=NULL)
{

  y <- n <- anc.values <- NULL
  m <-ncol(yy$xx) # number of traits
  trait_array<-array(data=NA, dim=(c(length(yy$xx[,1]), 4, m)))

  for (i in 1:m){
    trait_array[,1,i]<-yy$xx[,i]
    trait_array[,2,i]<-yy$vv[,i]
    trait_array[,3,i]<-yy$nn[,i]
    trait_array[,4,i]<-yy$tt[,i]
  }

  logic.init.diag.A<-diag(A.user)!=0
  nr.init.diag.A<-length(logic.init.diag.A[logic.init.diag.A==TRUE])
  init.diag.A<-rep(NA, nr.init.diag.A)
  locations.A<-which(A.user !=0, arr.ind = T)
  location.diag.A<-which(locations.A[,1] == locations.A[,2])

  location.upper.tri.A<-which(locations.A[,1] < locations.A[,2])
  nr.upper.tri.A<-length(location.upper.tri.A)

  location.lower.tri.A<-which(locations.A[,1] > locations.A[,2])
  nr.lower.tri.A<-length(location.lower.tri.A)

  nr.init.diag.R<-length(diag(R.user))
  init.diag.R<-rep(NA, nr.init.diag.R)

  locations.R<-which(R.user !=0, arr.ind = T)
  location.diag.R<-which(locations.R[,1] == locations.R[,2])

  location.upper.tri.R<-which(locations.R[,1] < locations.R[,2])
  nr.upper.tri.R<-length(location.upper.tri.R)

  for (i in 1:nr.init.diag.A)
  {
    init.diag.A[i]<-paleoTS::opt.joint.OU(paleoTS::as.paleoTS(mm=trait_array[,1,i], vv=trait_array[,2,i], nn=trait_array[,3,i], tt=trait_array[,4,i]))$parameter[4]
  }

  for (i in 1:nr.init.diag.R)
  {
    init.diag.R[i]<-paleoTS::opt.joint.URW(paleoTS::as.paleoTS(mm=trait_array[,1,i], vv=trait_array[,2,i], nn=trait_array[,3,i], tt=trait_array[,4,i]))$parameter[2]
  }

  init.upper.diag.A<-rep(-0.5, nr.upper.tri.A)
  init.lower.diag.A<-rep(0.5, nr.lower.tri.A)
  init.off.diag.R<-rep(0.5,nr.upper.tri.R)

  init.anc<-yy$xx[1,]

  init.theta<-init.anc
  tmp_diag_A.user<-diag(A.user)

  for (i in 1:m)
    {
    if (tmp_diag_A.user[i]==1) {init.theta[i]<-yy$xx[length(yy$xx[,1]),i]}
  }

  init.par<-c(init.diag.A, init.upper.diag.A, init.lower.diag.A, init.diag.R, init.off.diag.R, init.theta, init.anc)
  lower.limit<-c(rep(NA,length(init.diag.A)), rep(NA,length(init.upper.diag.A)),  rep(NA,length(init.lower.diag.A)), rep(0, length(init.diag.R)), rep(0, length(init.off.diag.R)), rep(NA, length(init.theta)), rep(NA, length(init.anc)))

  ##### Start of iteration routine #####

  if (is.numeric(iterations)) {
    if(is.numeric(iter.sd) == FALSE) iter.sd <-1
    log.lik.tmp<-rep(NA, iterations)
    www<-list()

    for (k in 1:iterations){

      #init.par<-rnorm(length(init.par_temp), init.par_temp, iter.sd)

  init.par_temp<-c(init.diag.A, init.upper.diag.A, init.lower.diag.A, init.diag.R, init.off.diag.R, init.theta, init.anc)
  init.par<-rnorm(length(init.par_temp), init.par_temp, iter.sd)
  lower.limit<-c(rep(NA,length(init.diag.A)), rep(NA,length(init.upper.diag.A)),  rep(NA,length(init.lower.diag.A)), rep(0, length(init.diag.R)), rep(0, length(init.off.diag.R)), rep(NA, length(init.theta)), rep(NA, length(init.anc)))


     if (method == "Nelder-Mead")  {
      www[[k]]<-optim(init.par, fn = logL.joint.multi.OUOU.user, yy = yy, A.user = A.user, R.user = R.user,
                      locations.A = locations.A, location.diag.A = location.diag.A, location.upper.tri.A = location.upper.tri.A, location.lower.tri.A = location.lower.tri.A,
                      locations.R = locations.R, location.diag.R = location.diag.R, location.upper.tri.R = location.upper.tri.R,
                       control = list(fnscale = -1, maxit=10000, trace = trace), method = "Nelder-Mead", hessian = hess)
     }

  if (method == "L-BFGS-B")  {
    www[[k]]<-optim(init.par, fn = logL.joint.multi.OUOU.user, yy = yy, A.user = A.user, R.user = R.user,
                    locations.A = locations.A, location.diag.A = location.diag.A, location.upper.tri.A = location.upper.tri.A, location.lower.tri.A = location.lower.tri.A,
                    locations.R = locations.R, location.diag.R = location.diag.R, location.upper.tri.R = location.upper.tri.R,
                    control = list(fnscale = -1, maxit=10000, trace = trace), method = "Nelder-Mead", hessian = hess, lower = lower.limit)
  }

    log.lik.tmp[k]<-www[[k]]$value

    }
    for (j in 1:iterations){
      if(max(na.exclude(log.lik.tmp)) == www[[j]]$value) best.run<-www[[j]]
    }
  }

  ##### End of iteration routine #####

  ##### Start of non-iteration routine #####

  if (is.numeric(iterations) == FALSE) {

    if (method == "Nelder-Mead")  {
      w<-optim(init.par, fn = logL.joint.multi.OUOU.user, yy = yy, A.user = A.user, R.user = R.user,
               locations.A = locations.A, location.diag.A = location.diag.A, location.upper.tri.A = location.upper.tri.A, location.lower.tri.A = location.lower.tri.A,
               locations.R = locations.R, location.diag.R = location.diag.R, location.upper.tri.R = location.upper.tri.R,
               control = list(fnscale = -1, maxit=10000, trace = trace), method = "Nelder-Mead", hessian = hess)
    }
    if (method == "L-BFGS-B")  {
      w<-optim(init.par, fn = logL.joint.multi.OUOU.user, yy = yy, A.user = A.user, R.user = R.user,
                locations.A = locations.A, location.diag.A = location.diag.A, location.upper.tri.A = location.upper.tri.A, location.lower.tri.A = location.lower.tri.A,
                locations.R = locations.R, location.diag.R = location.diag.R, location.upper.tri.R = location.upper.tri.R,
                control = list(fnscale = -1, maxit=10000, trace = trace), method = "Nelder-Mead", hessian = hess, lower = lower.limit)
    }

  }

  # number of parameters
  K <- nr.init.diag.A+nr.upper.tri.A+nr.lower.tri.A+nr.init.diag.R+nr.upper.tri.R+m+nr.init.diag.A

  if (is.numeric(iterations) == FALSE){
    iter<-NA
    if (hess) w$se <- sqrt(diag(-1 * solve(w$hessian))) else w$se <- NULL
  }

  if (is.numeric(iterations) == TRUE){
  iter<-iterations
  w<-best.run
  if (hess) best.run$se <- sqrt(diag(-1 * solve(w$hessian))) else w$se <- NULL
  }


  A<-diag(rep(0,m))
  A[locations.A[location.diag.A],locations.A[location.diag.A]]<- diag(c(w$par[1:length(location.diag.A)]))

  if (pracma::isempty(location.upper.tri.A)==FALSE)
  {
    A[locations.A[,1][location.upper.tri.A],locations.A[,2][location.upper.tri.A]]<-w$par[(length(location.diag.A)+1):(length(location.diag.A)+length(location.upper.tri.A))]
  } else location.upper.tri.A<-NULL

  if (pracma::isempty(location.lower.tri.A)==FALSE)
  {
    A[locations.A[,1][location.lower.tri.A],locations.A[,2][location.lower.tri.A]]<-w$par[(length(location.diag.A)+length(location.upper.tri.A)+1):(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A))]
  } else location.lower.tri.A<-NULL


  Chol<-diag(rep(0,m))
  Chol[locations.R[location.diag.R],locations.R[location.diag.R]]<- diag(c(w$par[(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+1):(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+length(location.diag.R))]))

  if (pracma::isempty(location.upper.tri.R)==FALSE)
  {
    Chol[locations.R[,1][location.upper.tri.R],locations.R[,2][location.upper.tri.R]]<-w$par[(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+length(location.diag.R)+1):(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+length(location.diag.R)+length(location.upper.tri.R))]
  } else location.upper.tri.R<-NULL
  R<-Chol %*% t(Chol)

  ### Theta (optimal trait values) ###
  optima<-c(w$par[(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+length(location.diag.R)+length(location.upper.tri.R)+1):(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+length(location.diag.R)+length(location.upper.tri.R)+m)])
#  if (nr.init.diag.A != 0) {
#    tmp_NA<-rep(NA, (m-nr.init.diag.A))
#    optima[(m-length(tmp_NA)+1):m]<-tmp_NA
#  }

  if (nr.init.diag.A < length(diag(A.user))) {
       tmp_NA<-rep(NA, (m-nr.init.diag.A))
      optima[(m-length(tmp_NA)+1):m]<-tmp_NA
      }

  ### The ancestral trait values ###
  ancestral.values<-c(w$par[(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+length(location.diag.R)+length(location.upper.tri.R)+m+1):(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+length(location.diag.R)+length(location.upper.tri.R)+m+m)])


  half.life<-log(2)/diag(A)


    wc<-as.evoTS.multi.OU.fit(logL = w$value, ancestral.values = ancestral.values, optima = optima, A = A, half.life = half.life, R = R,
                                                       method = "Joint", K = K, n = length(yy$xx[,1]), iter=iter)

  return(wc)

 }