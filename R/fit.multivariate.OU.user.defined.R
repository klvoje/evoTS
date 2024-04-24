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
#' @param pool indicating whether to pool variances across samples
#' 
#' @param trace logical, indicating whether information on the progress of the optimization is printed.
#'
#' @param iterations the number of times the optimization method is run from different starting points. Default is NULL, meaning the optimization is run once.
#'
#' @param iter.sd defines the standard deviation of the Gaussian distribution from which starting values for the optimization routine is run. Default is 1.
#'
#' @param user.init.diag.A starting values for the optimization routine of the diagonal elements of the A matrix. Default is NULL. 
#' 
#' @param user.init.upper.diag.A starting values for the optimization routine of the upper diagonal elements of the A matrix. Default is NULL.
#' 
#' @param user.init.lower.diag.A starting values for the optimization routine of the lower diagonal elements of the A matrix. Default is NULL.
#' 
#' @param user.init.diag.R starting values for the optimization routine of the diagonal elements of the R matrix. Default is NULL. 
#' 
#' @param user.init.off.diag.R starting values for the optimization routine of the off-diagonal elements of the R matrix. Default is NULL.
#' 
#' @param user.init.theta starting values for the optimization routine of the optima. Default is NULL.
#' 
#' @param user.init.anc starting values for the optimization routine of the ancestral values. Default is NULL.
#' 
#' @details This function provides users the flexibility to define their own A and R matrices. The possibility to define any A matrices enable detailed investigation of specific evolutionary hypotheses. The parameters to be estimated in the matrices are indicated by the value 1. All other entries in the matrix must be 0.
#'
#' The function searches - using an optimization routine - for the maximum-likelihood solution for the chosen multivariate Ornstein-Uhlenbeck model. The argument 'method' is passed to the 'optim' function and is included for the convenience of users to better control the optimization routine. Note that the the default method (Nelder-Mead) seems to work for most evolutionary sequences. The method L-BFGS-B allows box-constraints on some parameters (e.g. non-negative variance parameters) and is faster than Nelder-Mead, but is less stable than the default method (Nelder-Mead).
#'
#' Initial estimates to start the optimization come from maximum-likelihood estimates of the univariate Ornstein-Uhlenbeck model (from the paleoTS package) fitted to each time-series separately.
#'
#' It is good practice to repeat any numerical optimization procedure from different starting points. This is especially important for complex models as the log-likelihood surface might contain more than one peak. The number of iterations is controlled by the argument 'iterations'. The function will report the model parameters from the iteration with the highest log-likelihood.
#' 
#' There is no guarantee that the likelihood can be computed with the initial parameters provided by the function. The starting values for fitting the multivariate OU model are based on maximum likelihood parameter estimates for the univariate OU model fitted to each trait separately, which seems to provide sensible (and working) initial parameter estimates for almost all tested data sets. However, the provided initial parameters may fail depending on the nature of the data. If an error message is returned saying "function cannot be evaluated at initial parameters", the user can try to start the optimization procedure from other initial parameter values using "user.init.diag.A", "user.init.upper.diag.A", "user.init.lower.diag.A", "user.init.diag.R", "user.init.off.diag.R", "user.init.theta", and "user.init.anc." It is usually the initial guess of the off-diagonal elements of the A and R matrices that prevents the optimization routine to work. It is therefore recommended to only try to change these initial values before experimenting with different starting values for the diagonal of the A and R matrices.  
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
#'
#'## Generate a evoTS object by simulating a multivariate dataset
#'x <- sim.multi.OU(15)
#'
#'## Define an A matrix that is lower diagonal.
#'A <- matrix(c(1,0,1,1), nrow=2, byrow=TRUE)
#'
#'## Define a diagonal R matrix.
#'R <- matrix(c(1,0,0,1), nrow=2, byrow=TRUE)
#'
#'\donttest{
#'## Fit the multivariate Ornstein-Uhlenbeck model to the data. This example will run for a long time.
#'fit.multivariate.OU.user.defined(x, A.user=A, R.user=R, trace=TRUE)
#'}

fit.multivariate.OU.user.defined<-function (yy, A.user=NULL, R.user=NULL, method="Nelder-Mead", hess = FALSE, pool = TRUE, trace=FALSE, iterations=NULL, iter.sd=NULL, user.init.diag.A = NULL, user.init.upper.diag.A = NULL, user.init.lower.diag.A = NULL, user.init.diag.R = NULL, user.init.off.diag.R = NULL, user.init.theta = NULL, user.init.anc = NULL)
  {

  y <- n <- anc.values <- NULL
  m <-ncol(yy$xx) # number of traits
  
  
  if (pool==TRUE) { 
    for (i in 1:m){
      
      tmp<-paleoTS::as.paleoTS(yy$xx[,i], yy$vv[,i], yy$nn[,i], yy$tt[,i])
      tmp<- paleoTS::pool.var(tmp, ret.paleoTS = TRUE)
      yy$vv[,i]<-tmp$vv
    }
  }
  
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

  init.upper.diag.A<-rep(0, nr.upper.tri.A)
  init.lower.diag.A<-rep(0, nr.lower.tri.A)
  init.off.diag.R<-rep(0.5,nr.upper.tri.R)

  init.anc<-yy$xx[1,]

  init.theta<-init.anc
  tmp_diag_A.user<-diag(A.user)

  for (i in 1:m)
    {
    if (tmp_diag_A.user[i]==1) {init.theta[i]<-yy$xx[length(yy$xx[,1]),i]}
  }


  #Check for user defined starting values
  if (is.null(user.init.diag.A) == FALSE) init.diag.A<-user.init.diag.A
  if (is.null(user.init.upper.diag.A) == FALSE) init.upper.diag.A<-user.init.upper.diag.A
  if (is.null(user.init.lower.diag.A) == FALSE) init.lower.diag.A<-user.init.lower.diag.A
  if (is.null(user.init.diag.R) == FALSE) init.diag.R<-user.init.diag.R
  if (is.null(user.init.off.diag.R) == FALSE) init.off.diag.R<-user.init.off.diag.R
  if (is.null(user.init.theta) == FALSE) init.theta<-user.init.theta
  if (is.null(user.init.anc) == FALSE) init.anc<-user.init.anc

  init.par<-c(init.diag.A, init.upper.diag.A, init.lower.diag.A, init.diag.R, init.off.diag.R, init.theta, init.anc)
  lower.limit<-c(rep(NA,length(init.diag.A)), rep(NA,length(init.upper.diag.A)),  rep(NA,length(init.lower.diag.A)), rep(0, length(init.diag.R)), rep(0, length(init.off.diag.R)), rep(NA, length(init.theta)), rep(NA, length(init.anc)))

  ##### Start of iteration routine #####

  if (is.numeric(iterations)) {
    if(is.numeric(iter.sd) == FALSE) iter.sd <-1
    log.lik.tmp<-rep(NA, 1000000)
    www<-list()

    for (k in 1:1000000){
      tryCatch({
      #init.par<-rnorm(length(init.par_temp), init.par_temp, iter.sd)

  init.par_temp<-c(init.diag.A, init.upper.diag.A, init.lower.diag.A, init.diag.R, init.off.diag.R, init.theta, init.anc)
  init.par<-rnorm(length(init.par_temp), init.par_temp, iter.sd)
  lower.limit<-c(rep(NA,length(init.diag.A)), rep(NA,length(init.upper.diag.A)),  rep(NA,length(init.lower.diag.A)), rep(0, length(init.diag.R)), rep(0, length(init.off.diag.R)), rep(NA, length(init.theta)), rep(NA, length(init.anc)))


     if (method == "Nelder-Mead")  {
      www[[k]]<-try(optim(init.par, fn = logL.joint.multi.OUOU.user, yy = yy, A.user = A.user, R.user = R.user,
                      locations.A = locations.A, location.diag.A = location.diag.A, location.upper.tri.A = location.upper.tri.A, location.lower.tri.A = location.lower.tri.A,
                      locations.R = locations.R, location.diag.R = location.diag.R, location.upper.tri.R = location.upper.tri.R,
                       control = list(fnscale = -1, maxit=1000000, trace = trace), method = "Nelder-Mead", hessian = hess), silent = TRUE)
      if(inherits(www[[k]], "try-error") && grepl("function cannot be evaluated at initial parameters", attr(www[[k]], "condition")$message))
           stop("The initial parameters did not work. Trying a new set of candidate starting values.")
      # The provided initial starting values for the parameters may not work (depends on the data). If this happens when running iterations, the user is informed by a message saying: "The initial parameters did not work. Trying a new set of candidate starting values." 
      
     }
    
  if (method == "L-BFGS-B")  {
    www[[k]]<-optim(init.par, fn = logL.joint.multi.OUOU.user, yy = yy, A.user = A.user, R.user = R.user,
                    locations.A = locations.A, location.diag.A = location.diag.A, location.upper.tri.A = location.upper.tri.A, location.lower.tri.A = location.lower.tri.A,
                    locations.R = locations.R, location.diag.R = location.diag.R, location.upper.tri.R = location.upper.tri.R,
                    control = list(fnscale = -1, maxit=1000000, trace = trace), method = "L-BFGS-B", hessian = hess, lower = lower.limit)
  }
    
    log.lik.tmp[k]<-www[[k]]$value

    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

        if (length(na.exclude(log.lik.tmp)) == iterations){
        break
      }
    }
    
    # Need to remove entries in www in case there are iterations where the initial parameter estimates did not work.
    www_tmp<-list()
    for (i in 1:k){
      if (is.character(www[[i]][1]) == FALSE) www_tmp[i]<-list(www[[i]])
    }
    
    www_reduced<-www_tmp[!sapply(www_tmp,is.null)]
    
    
    for (j in 1:length(www_reduced)){
      if(max(na.exclude(log.lik.tmp)) == www_reduced[[j]]$value) best.run<-www_reduced[[j]]
    }
    
    }

  
  ##### End of iteration routine #####

  ##### Start of non-iteration routine #####

  if (is.numeric(iterations) == FALSE) {

    if (method == "Nelder-Mead")  {
      w<-optim(init.par, fn = logL.joint.multi.OUOU.user, yy = yy, A.user = A.user, R.user = R.user,
               locations.A = locations.A, location.diag.A = location.diag.A, location.upper.tri.A = location.upper.tri.A, location.lower.tri.A = location.lower.tri.A,
               locations.R = locations.R, location.diag.R = location.diag.R, location.upper.tri.R = location.upper.tri.R,
               control = list(fnscale = -1, maxit=1000000, trace = trace), method = "Nelder-Mead", hessian = hess)
    }
    if (method == "L-BFGS-B")  {
      w<-optim(init.par, fn = logL.joint.multi.OUOU.user, yy = yy, A.user = A.user, R.user = R.user,
                locations.A = locations.A, location.diag.A = location.diag.A, location.upper.tri.A = location.upper.tri.A, location.lower.tri.A = location.lower.tri.A,
                locations.R = locations.R, location.diag.R = location.diag.R, location.upper.tri.R = location.upper.tri.R,
                control = list(fnscale = -1, maxit=1000000, trace = trace), method = "Nelder-Mead", hessian = hess, lower = lower.limit)
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
  if (hess) w$se <- sqrt(diag(-1 * solve(w$hessian))) else w$se <- NULL
  }

  if (w$convergence == 10) converge<-"The search algorithm stopped as it did not make progress towards the optimal solution"
  if (w$convergence == 0) converge<-"Model converged successfully"
  if (w$convergence == 1) converge<-"The maximum number of iterations was reached and the search algorithm exited"
  
  A<-diag(rep(0,m))

  for (i in 1:length(location.diag.A)){
    A[locations.A[location.diag.A][i],locations.A[location.diag.A][i]]<- w$par[i]
  }

  if (pracma::isempty(location.upper.tri.A)==FALSE)
  {
    
    for (i in 1:length(location.upper.tri.A)){
      A[locations.A[,1][location.upper.tri.A][i],locations.A[,2][location.upper.tri.A][i]]<- w$par[(length(location.diag.A)+i)]
    }
  } else location.upper.tri.A<-NULL
  

  if (pracma::isempty(location.lower.tri.A)==FALSE)
  {
    for (i in 1:length(location.lower.tri.A)){
      A[locations.A[,1][location.lower.tri.A][i],locations.A[,2][location.lower.tri.A][i]]<-w$par[(length(location.diag.A)+length(location.upper.tri.A)+i)]
    } 
  }else location.lower.tri.A<-NULL


  Chol<-diag(rep(0,m))
  for (i in 1:length(location.diag.R)){
    Chol[locations.R[location.diag.R][i],locations.R[location.diag.R][i]]<- w$par[(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+i)]
  }
  
  if (pracma::isempty(location.upper.tri.R)==FALSE)
  {
    for (i in 1:length(location.upper.tri.R)){
      Chol[locations.R[,1][location.upper.tri.R][i],locations.R[,2][location.upper.tri.R][i]]<-w$par[(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+length(location.diag.R)+i)]
    }
  } else location.upper.tri.R<-NULL
  
  R<-Chol %*% t(Chol)

   ### Theta (optimal trait values) ###
  optima<-c(w$par[(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+length(location.diag.R)+length(location.upper.tri.R)+1):(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+length(location.diag.R)+length(location.upper.tri.R)+m)])

  for (i in 1:m)
  {
    if (init.theta[i]==init.anc[i]) optima[i]<-NA 
  }


  ### The ancestral trait values ###
  ancestral.values<-c(w$par[(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+length(location.diag.R)+length(location.upper.tri.R)+m+1):(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+length(location.diag.R)+length(location.upper.tri.R)+m+m)])

  if (hess == TRUE){
  # SE parameters
  
  SE.A<-diag(rep(0,m))
  
  for (i in 1:length(location.diag.A)){
    SE.A[locations.A[location.diag.A][i],locations.A[location.diag.A][i]]<- w$se[i]
  }
  
  if (pracma::isempty(location.upper.tri.A)==FALSE)
  {
    
    for (i in 1:length(location.upper.tri.A)){
      SE.A[locations.A[,1][location.upper.tri.A][i],locations.A[,2][location.upper.tri.A][i]]<- w$se[(length(location.diag.A)+i)]
    }
  } else location.upper.tri.A<-NULL
  
  
  if (pracma::isempty(location.lower.tri.A)==FALSE)
  {
    for (i in 1:length(location.lower.tri.A)){
      SE.A[locations.A[,1][location.lower.tri.A][i],locations.A[,2][location.lower.tri.A][i]]<-w$se[(length(location.diag.A)+length(location.upper.tri.A)+i)]
    } 
  }else location.lower.tri.A<-NULL
  
  
  SE.Chol<-diag(rep(0,m))
  for (i in 1:length(location.diag.R)){
    SE.Chol[locations.R[location.diag.R][i],locations.R[location.diag.R][i]]<- w$se[(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+i)]
  }
  
  if (pracma::isempty(location.upper.tri.R)==FALSE)
  {
    for (i in 1:length(location.upper.tri.R)){
      SE.Chol[locations.R[,1][location.upper.tri.R][i],locations.R[,2][location.upper.tri.R][i]]<-w$se[(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+length(location.diag.R)+i)]
    }
  } else location.upper.tri.R<-NULL
  
  SE.R<-SE.Chol %*% t(SE.Chol)
  
  ### Theta (optimal trait values) ###
  SE.optima<-c(w$se[(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+length(location.diag.R)+length(location.upper.tri.R)+1):(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+length(location.diag.R)+length(location.upper.tri.R)+m)])
  
  for (i in 1:m)
  {
    if (init.theta[i]==init.anc[i]) SE.optima[i]<-NA 
  }
  
  ### The ancestral trait values ###
  SE.anc<-c(w$se[(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+length(location.diag.R)+length(location.upper.tri.R)+m+1):(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+length(location.diag.R)+length(location.upper.tri.R)+m+m)])
  }
  
  half.life<-log(2)/diag(A)
  
  modelName<-"Multivariate model: User-specified OU model"

  if (hess == TRUE){
    
    wc<-as.evoTS.multi.OU.fit(converge, modelName = modelName, logL = w$value, ancestral.values = ancestral.values, SE.anc = SE.anc, optima = optima,  SE.optima = SE.optima, A = A, SE.A = SE.A, half.life = half.life, R = R, SE.R = SE.R,
                              method = "Joint", K = K, n = length(yy$xx[,1]), iter=iter)
  } 
  
  if (hess == FALSE){
    SE.anc <- NA
    SE.optima <- NA
    SE.A <- NA 
    SE.R <- NA 
    wc<-as.evoTS.multi.OU.fit(converge, modelName = modelName, logL = w$value, ancestral.values = ancestral.values, SE.anc = SE.anc, optima = optima,  SE.optima = SE.optima, A = A, SE.A = SE.A, half.life = half.life, R = R, SE.R = SE.R,
                                                       method = "Joint", K = K, n = length(yy$xx[,1]), iter=iter)
  }
  
  return(wc)

 }
