#' @title Fit predefined multivariate Ornstein-Uhlenbeck models to multivariate evolutionary sequence (time-series) data.
#'
#' @description Function to find maximum likelihood solutions to a large suite of predefined multivariate Ornstein-Uhlenbeck model fitted to multivariate evolutionary sequence (time-series) data.
#'
#' @param yy a multivariate evoTS object.
#'
#' @param A.matrix the pull matrix. The options are "diag", "upper.tri", "lower.tri", and "full". Default is "diag". See details (or vignette) for more info on what the different options mean.
#'
#' @param R.matrix the drift matrix. The options are "diag" and "symmetric".
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
#' @param user.init.diag.R starting values for the optimization routine of the diagonal elements of the R matrix. Default is NULL.
#' 
#' @param user.init.off.diag.A starting values for the optimization routine of the off-diagonal elements of the A matrix. Default is NULL.
#' 
#' @param user.init.off.diag.R starting values for the optimization routine of the off-diagonal elements of the R matrix. Default is NULL.
#' 
#' @param user.init.theta starting values for the optimization routine of the optima. Default is NULL.
#' 
#' @param user.init.anc starting values for the optimization routine of the ancestral values. Default is NULL.
#'
#' @details A detailed explanation of the predefined models that can be fitted using the function is given in the online vignette (https://klvoje.github.io/evoTS/index.html), but a short summary is provided here. Note that this function provides the user with fixed options for how to parameterize the A and R matrices. For full flexibility, the user is allowed to customize the parameterization of the A and R matrix in the 'fit.multivariate.OU.user.defined' function.
#' The type of trait dynamics is defined based on how the pull matrix (A) and drift matrix (R) are defined. The function allows testing four broad categories of models: 1 Independent evolution (A.matrix ="diag", R.matrix = "diag"); 2 Independent adaptation (A.matrix ="diag", R.matrix = "symmetric"); 3 Non-independent adaptation (A.matrix = "upper.tri"/"lower.tri"/full", R.matrix = "diagonal"); 4 Non-independent evolution (A.matrix = "upper.tri"/"lower.tri"/"full", R.matrix = "symmetric").
#' Setting the A.matrix to "diagonal" means the traits do not affect each others optimum (A matrix). A "diagonal" R matrix means the stochastic changes in the traits are assumed to be uncorrelated. A "symmetric" R matrix means the stochastic changes in the traits are assumed to be correlated, i.e. that they are non-independent. A "full" parameterization of A estimates the effect of each trait on the optima on the other traits.
#' The "upper.tri" option parameterize the model in such a way that the first layer (first trait in the data set) adapts non-independently because its optimum is affected by all other traits included in the data set, while the bottom layer (the last trait in the data set) adapts independently (as an Ornstein Uhlenbeck process). Layers in between the upper- and lower layer (not the first or last trait in the data set (if there are more than two traits in the data set)) evolve non-independently as their optimum is affected by all layers/traits below themselves. The option "lower.tri" defines the causality the opposite way compared to "upper.tri".
#' It is also possible to implement a model where the bottom layer (last trait in the data set) evolve as an Unbiased random walk (akin to a Brownian motion) which affects the optima for all other traits in the data set (i.e. all layers except the bottom layer). This model can be fitted by defining A.matrix ="OUBM", which will override how the R matrix is defined.
#'
#' The function searches - using an optimization routine - for the maximum-likelihood solution for the chosen multivariate Ornstein-Uhlenbeck model. The argument 'method' is passed to the 'optim' function and is included for the convenience of users to better control the optimization routine. Note that the the default method (Nelder-Mead) seems to work for most evolutionary sequences. The method L-BFGS-B allows box-constraints on some parameters (e.g. non-negative variance parameters) and is faster than Nelder-Mead, but is less stable than the default method (Nelder-Mead).
#'
#' Initial estimates to start the optimization come from maximum-likelihood estimates of the univariate Ornstein-Uhlenbeck model (from the paleoTS package) fitted to each time-series separately.
#'
#' It is good practice to repeat any numerical optimization procedure from different starting points. This is especially important for complex models as the log-likelihood surface might contain more than one peak. The number of iterations is controlled by the argument 'iterations'. The function will report the model parameters from the iteration with the highest log-likelihood.
#' 
#' There is no guarantee that the likelihood can be computed with the initial parameters provided by the function. The starting values for fitting the multivariate OU model are based on maximum likelihood parameter estimates for the univariate OU model fitted to each trait separately, which seems to provide sensible (and working) initial parameter estimates for almost all tested data sets. However, the provided initial parameters may fail depending on the nature of the data. If an error message is returned saying "function cannot be evaluated at initial parameters", the user can try to start the optimization procedure from other initial parameter values using "user.init.diag.A", "user.init.diag.R", "user.init.off.diag.A", "user.init.off.diag.R", "user.init.theta", and "user.init.anc." It is usually the initial guess of the off-diagonal elements of the A and R matrices that prevents the optimization routine to work. It is therefore recommended to only try to change these initial values before experimenting with different starting values for the diagonal of the A and R matrices.  
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
#'\donttest{
#'##Fit a multivariate Ornstein-Uhlenbeck model to the data. This example will run for a long time.
#'fit.multivariate.OU(x, A.matrix="diag", R.matrix="symmetric")
#'}

fit.multivariate.OU<-function (yy, A.matrix="diag", R.matrix="symmetric", method="Nelder-Mead", hess = FALSE, pool = TRUE, trace=FALSE, iterations=NULL, iter.sd=NULL, user.init.diag.A = NULL, user.init.diag.R = NULL, user.init.off.diag.A = NULL, user.init.off.diag.R = NULL, user.init.theta = NULL, user.init.anc = NULL)
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
    
  
  ### Start iterations from different starting values
  if (is.numeric(iterations)) {
    if(is.numeric(iter.sd) == FALSE) iter.sd <-1
    log.lik.tmp<-rep(NA, 1000000)
    www<-list()

    for (k in 1:1000000){
      tryCatch({
      trait_array<-array(data=NA, dim=(c(length(yy$xx[,1]), 4, m)))

      for (i in 1:m){
        trait_array[,1,i]<-yy$xx[,i]
        trait_array[,2,i]<-yy$vv[,i]
        trait_array[,3,i]<-yy$nn[,i]
        trait_array[,4,i]<-yy$tt[,i]
      }

      init.diag.A<-rep(NA,m)
      init.diag.R<-init.diag.A

      if (A.matrix=="OUBM"){init.diag.A<-rep(NA,(m-1))}

      for (i in 1:length(init.diag.A))
      {
        init.diag.A[i]<-paleoTS::opt.joint.OU(paleoTS::as.paleoTS(mm=trait_array[,1,i], vv=trait_array[,2,i], nn=trait_array[,3,i], tt=trait_array[,4,i]))$parameter[4]
      }

      for (i in 1:length(init.diag.R))
      {
        init.diag.R[i]<-paleoTS::opt.joint.URW(paleoTS::as.paleoTS(mm=trait_array[,1,i], vv=trait_array[,2,i], nn=trait_array[,3,i], tt=trait_array[,4,i]))$parameter[2]
      }


      init.off.diag.A<-rep(0, sum(upper.tri(diag(init.diag.A)), na.rm = TRUE))
      if (A.matrix=="OUBM"){
        if (length(init.diag.A) == 1) init.off.diag.A<-0
        if (length(init.diag.A) != 1) init.off.diag.A<-rep(0, sum(upper.tri(diag(init.diag.A)), na.rm = TRUE))
        }
      if (A.matrix=="full"){init.off.diag.A<-rep(0, (sum(upper.tri(diag(init.diag.A)), na.rm = TRUE))*2)}
      init.off.diag.R<-rep(0.5, sum(upper.tri(diag(init.diag.R)), na.rm = TRUE))

      init.anc<-yy$xx[1,]

      init.theta<-yy$xx[length(yy$xx[,1]),]
      if (A.matrix=="OUBM"){init.theta[m]<-init.anc[m]}
      
      #Check for user defined starting values
      if (is.null(user.init.diag.A) == FALSE) init.diag.A<-user.init.diag.A
      if (is.null(user.init.diag.R) == FALSE) init.diag.R<-user.init.diag.R
      if (is.null(user.init.off.diag.A) == FALSE) init.off.diag.A<-user.init.off.diag.A
      if (is.null(user.init.off.diag.R) == FALSE) init.off.diag.R<-user.init.off.diag.R
      if (is.null(user.init.theta) == FALSE) init.theta<-user.init.theta
      if (is.null(user.init.anc) == FALSE) init.anc<-user.init.anc
      
    if(A.matrix=="diag" & R.matrix=="diag")
    {
      init.par_temp<-c(init.diag.A, init.diag.R, init.theta, init.anc)
      init.par<-rnorm(length(init.par_temp), init.par_temp, iter.sd)
      lower.limit<-c(rep(NA,length(init.diag.A)), rep(0, length(init.diag.R)), rep(NA, length(init.theta)), rep(NA, length(init.anc)))
    }

    if(A.matrix=="diag" & R.matrix=="symmetric")
    {
      init.par_temp<-c(init.diag.A, init.diag.R, init.off.diag.R, init.theta, init.anc)
      init.par<-rnorm(length(init.par_temp), init.par_temp, iter.sd)
      lower.limit<-c(rep(NA,length(init.diag.A)), rep(0, length(init.diag.R)), rep(NA, length(init.off.diag.R)), rep(NA, length(init.theta)), rep(NA, length(init.anc)))
    }

      if(A.matrix=="full" & R.matrix=="symmetric")
      {
        init.par_temp<-c(init.diag.A, init.off.diag.A, init.diag.R, init.off.diag.R,  init.theta, init.anc)
        init.par<-rnorm(length(init.par_temp), init.par_temp, iter.sd)
        lower.limit<-c(rep(NA,length(init.diag.A)), rep(NA,length(init.off.diag.A)), rep(0, length(init.diag.R)), rep(NA, length(init.off.diag.R)), rep(NA, length(init.theta)), rep(NA, length(init.anc)))
      }

      if(A.matrix=="full" & R.matrix=="diag")
      {
        init.par_temp<-c(init.diag.A, init.off.diag.A, init.diag.R,  init.theta, init.anc)
        init.par<-rnorm(length(init.par_temp), init.par_temp, iter.sd)
        lower.limit<-c(rep(NA,length(init.diag.A)), rep(NA,length(init.off.diag.A)), rep(0, length(init.diag.R)), rep(NA, length(init.theta)), rep(NA, length(init.anc)))
      }

    if(A.matrix=="upper.tri" & R.matrix=="diag")
    {
      init.par_temp<-c(init.diag.A, init.off.diag.A, init.diag.R,  init.theta, init.anc)
      init.par<-rnorm(length(init.par_temp), init.par_temp, iter.sd)
      lower.limit<-c(rep(NA,length(init.diag.A)), rep(NA,length(init.off.diag.A)), rep(0, length(init.diag.R)), rep(NA, length(init.theta)), rep(NA, length(init.anc)))
    }

    if(A.matrix=="upper.tri" & R.matrix=="symmetric")
    {
      init.par_temp<-c(init.diag.A, init.off.diag.A, init.diag.R, init.off.diag.R,  init.theta, init.anc)
      init.par<-rnorm(length(init.par_temp), init.par_temp, iter.sd)
      lower.limit<-c(rep(NA,length(init.diag.A)), rep(NA,length(init.off.diag.A)), rep(0, length(init.diag.R)), rep(NA, length(init.off.diag.R)), rep(NA, length(init.theta)), rep(NA, length(init.anc)))
    }

      if(A.matrix=="lower.tri" & R.matrix=="diag")
      {
        init.par_temp<-c(init.diag.A, init.off.diag.A, init.diag.R,  init.theta, init.anc)
        init.par<-rnorm(length(init.par_temp), init.par_temp, iter.sd)
        lower.limit<-c(rep(NA,length(init.diag.A)), rep(NA,length(init.off.diag.A)), rep(0, length(init.diag.R)), rep(NA, length(init.theta)), rep(NA, length(init.anc)))
      }

      if(A.matrix=="lower.tri" & R.matrix=="symmetric")
      {
        init.par_temp<-c(init.diag.A, init.off.diag.A, init.diag.R, init.off.diag.R,  init.theta, init.anc)
        init.par<-rnorm(length(init.par_temp), init.par_temp, iter.sd)
        lower.limit<-c(rep(NA,length(init.diag.A)), rep(NA,length(init.off.diag.A)), rep(0, length(init.diag.R)), rep(NA, length(init.off.diag.R)), rep(NA, length(init.theta)), rep(NA, length(init.anc)))
      }

    if(A.matrix=="OUBM")
    {
      init.par_temp<-c(init.diag.A, init.off.diag.A, init.diag.R, init.off.diag.R,  init.theta, init.anc)
      init.par<-rnorm(length(init.par_temp), init.par_temp, iter.sd)
      lower.limit<-c(rep(NA,length(init.diag.A)), rep(NA,length(init.off.diag.A)), rep(0, length(init.diag.R)), rep(NA, length(init.off.diag.R)), rep(NA, length(init.theta)), rep(NA, length(init.anc)))
    }


     if (method == "Nelder-Mead")  {
      www[[k]]<-try(optim(init.par, fn = logL.joint.multi.OUOU, yy = yy, A.matrix = A.matrix, R.matrix = R.matrix,
                       control = list(fnscale = -1, maxit=1000000, trace = trace), method = "Nelder-Mead", hessian = hess), silent = TRUE)
      if(inherits(www[[k]], "try-error") && grepl("function cannot be evaluated at initial parameters", attr(www[[k]], "condition")$message))
        stop("The initial parameters did not work. Trying a new set of candidate starting values.")
      # The provided initial starting values for the parameters may not work (depends on the data). If this happens when running iterations, the user is informed by a message saying: "The initial parameters did not work. Trying a new set of candidate starting values." 
     }
      
      
    if (method == "L-BFGS-B")  {
      www[[k]]<-optim(init.par, fn = logL.joint.multi.OUOU, yy = yy, A.matrix = A.matrix, R.matrix = R.matrix,
                 control = list(fnscale = -1, maxit=1000000, trace = trace), method = "L-BFGS-B" , hessian = hess, lower = lower.limit)
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


  ##### Start of no-iteration routine #####

  if (is.numeric(iterations) == FALSE) {
    trait_array<-array(data=NA, dim=(c(length(yy$xx[,1]), 4, m)))

  for (i in 1:m){
    trait_array[,1,i]<-yy$xx[,i]
    trait_array[,2,i]<-yy$vv[,i]
    trait_array[,3,i]<-yy$nn[,i]
    trait_array[,4,i]<-yy$tt[,i]
  }

      init.diag.A<-rep(NA,m)
      init.diag.R<-init.diag.A

      if (A.matrix=="OUBM"){init.diag.A<-rep(NA,(m-1))}

        for (i in 1:length(init.diag.A))
          {
          init.diag.A[i]<-paleoTS::opt.joint.OU(paleoTS::as.paleoTS(mm=trait_array[,1,i], vv=trait_array[,2,i], nn=trait_array[,3,i], tt=trait_array[,4,i]))$parameter[4]
        }

         for (i in 1:length(init.diag.R))
          {
          init.diag.R[i]<-paleoTS::opt.joint.URW(paleoTS::as.paleoTS(mm=trait_array[,1,i], vv=trait_array[,2,i], nn=trait_array[,3,i], tt=trait_array[,4,i]))$parameter[2]
        }


      init.off.diag.A<-rep(0, sum(upper.tri(diag(init.diag.A)), na.rm = TRUE))
      if (A.matrix=="OUBM"){
        if (length(init.diag.A) == 1) init.off.diag.A<-0
        if (length(init.diag.A) != 1) init.off.diag.A<-rep(0, sum(upper.tri(diag(init.diag.A)), na.rm = TRUE))
      }
      if (A.matrix=="full"){init.off.diag.A<-rep(0, (sum(upper.tri(diag(init.diag.A)), na.rm = TRUE))*2)}
      init.off.diag.R<-rep(0.5, sum(upper.tri(diag(init.diag.R)), na.rm = TRUE))

      init.anc<-yy$xx[1,]

      init.theta<-yy$xx[length(yy$xx[,1]),]
      if (A.matrix=="OUBM"){init.theta[m]<-init.anc[m]}
      
      #Check for user defined starting values
      if (is.null(user.init.diag.A) == FALSE) init.diag.A<-user.init.diag.A
      if (is.null(user.init.diag.R) == FALSE) init.diag.R<-user.init.diag.R
      if (is.null(user.init.off.diag.A) == FALSE) init.off.diag.A<-user.init.off.diag.A
      if (is.null(user.init.off.diag.R) == FALSE) init.off.diag.R<-user.init.off.diag.R
      if (is.null(user.init.theta) == FALSE) init.theta<-user.init.theta
      if (is.null(user.init.anc) == FALSE) init.anc<-user.init.anc


      if(A.matrix=="diag" & R.matrix=="diag")
        {
      init.par<-c(init.diag.A, init.diag.R, init.theta, init.anc)
      lower.limit<-c(rep(NA,length(init.diag.A)), rep(0, length(init.diag.R)), rep(NA, length(init.theta)), rep(NA, length(init.anc)))
      }

      if(A.matrix=="diag" & R.matrix=="symmetric")
        {
        init.par<-c(init.diag.A, init.diag.R, init.off.diag.R, init.theta, init.anc)
        lower.limit<-c(rep(NA,length(init.diag.A)), rep(0, length(init.diag.R)), rep(NA, length(init.off.diag.R)), rep(NA, length(init.theta)), rep(NA, length(init.anc)))
      }

      if(A.matrix=="full" & R.matrix=="symmetric")
      {
        init.par<-c(init.diag.A, init.off.diag.A, init.diag.R, init.off.diag.R,  init.theta, init.anc)
        lower.limit<-c(rep(NA,length(init.diag.A)), rep(NA,length(init.off.diag.A)), rep(0, length(init.diag.R)), rep(NA, length(init.off.diag.R)), rep(NA, length(init.theta)), rep(NA, length(init.anc)))
      }

      if(A.matrix=="full" & R.matrix=="diag")
      {
        init.par<-c(init.diag.A, init.off.diag.A, init.diag.R,  init.theta, init.anc)
        lower.limit<-c(rep(NA,length(init.diag.A)), rep(NA,length(init.off.diag.A)), rep(0, length(init.diag.R)), rep(NA, length(init.theta)), rep(NA, length(init.anc)))
      }

      if(A.matrix=="upper.tri" & R.matrix=="diag")
      {
        init.par<-c(init.diag.A, init.off.diag.A, init.diag.R,  init.theta, init.anc)
        lower.limit<-c(rep(NA,length(init.diag.A)), rep(NA,length(init.off.diag.A)), rep(0, length(init.diag.R)), rep(NA, length(init.theta)), rep(NA, length(init.anc)))
      }

      if(A.matrix=="upper.tri" & R.matrix=="symmetric")
      {
        init.par<-c(init.diag.A, init.off.diag.A, init.diag.R, init.off.diag.R,  init.theta, init.anc)
        lower.limit<-c(rep(NA,length(init.diag.A)), rep(NA,length(init.off.diag.A)), rep(0, length(init.diag.R)), rep(NA, length(init.off.diag.R)), rep(NA, length(init.theta)), rep(NA, length(init.anc)))
      }

      if(A.matrix=="lower.tri" & R.matrix=="diag")
      {
        init.par<-c(init.diag.A, init.off.diag.A, init.diag.R,  init.theta, init.anc)
        lower.limit<-c(rep(NA,length(init.diag.A)), rep(NA,length(init.off.diag.A)), rep(0, length(init.diag.R)), rep(NA, length(init.theta)), rep(NA, length(init.anc)))
      }

      if(A.matrix=="lower.tri" & R.matrix=="symmetric")
      {
        init.par<-c(init.diag.A, init.off.diag.A, init.diag.R, init.off.diag.R,  init.theta, init.anc)
        lower.limit<-c(rep(NA,length(init.diag.A)), rep(NA,length(init.off.diag.A)), rep(0, length(init.diag.R)), rep(NA, length(init.off.diag.R)), rep(NA, length(init.theta)), rep(NA, length(init.anc)))
      }

      if(A.matrix=="OUBM")
      {
        init.par<-c(init.diag.A, init.off.diag.A, init.diag.R,  init.theta, init.anc)
        lower.limit<-c(rep(NA,length(init.diag.A)), rep(NA,length(init.off.diag.A)), rep(0, length(init.diag.R)), rep(NA, length(init.off.diag.R)), rep(NA, length(init.theta)), rep(NA, length(init.anc)))
      }


    if (method == "L-BFGS-B")
    {
      w<-optim(init.par, fn = logL.joint.multi.OUOU, yy = yy, A.matrix = A.matrix, R.matrix = R.matrix,
               control = list(fnscale = -1, maxit=1000000, trace = trace), method = "L-BFGS-B", hessian = hess, lower = lower.limit)
    }

    if (method == "Nelder-Mead")
    {
      w<-optim(init.par, fn = logL.joint.multi.OUOU, yy = yy, A.matrix = A.matrix, R.matrix = R.matrix,
               control = list(fnscale = -1, maxit=1000000, trace = trace), method = "Nelder-Mead", hessian = hess)
    }
  }

    ##### End of no-iteration routine #####


  # number of parameters
  K <- length(init.par)

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

  
  if(A.matrix=="diag" & R.matrix=="diag") {
    A<-diag(c(w$par[1:m]))

    Chol<-diag(c(w$par[(m+1):(m*2)]))
    R<-t(Chol)%*%Chol

    optima<-c(w$par[((m*2)+1):(m*3)])

    ancestral.values<-c(w$par[((m*3)+1):(m*4)])
    
    # SE parameters
    if (hess){
    SE.A<-diag(c(w$se[1:m]))
    
    SE.Chol<-diag(c(w$se[(m+1):(m*2)]))
    SE.R<-t(SE.Chol)%*%SE.Chol
    
    SE.optima<-c(w$se[((m*2)+1):(m*3)])
   
    SE.anc<-c(w$se[((m*3)+1):(m*4)])
    }
    
    modelName<-"Multivariate model: OU; diagonal A matrix, diagonal R matrix"
    
  }

  if(A.matrix=="diag" & R.matrix=="symmetric") {
    A<-diag(c(w$par[1:m]))
    Chol<-diag(w$par[(m+1):(m*2)])
    nr.off.diag<-upper.tri(Chol)
    l.upp.tri<-length(nr.off.diag[nr.off.diag==TRUE])

    Chol[upper.tri(Chol)] <- w$par[((length(diag(A))*2)+1):((length(diag(A))*2)+l.upp.tri)]
    R<-t(Chol)%*%Chol

    optima<-c(w$par[((length(diag(A))*2)+l.upp.tri+1):((length(diag(A))*2)+l.upp.tri+m)])

    ancestral.values<-c(w$par[((length(diag(A))*2)+l.upp.tri+1+m):((length(diag(A))*2)+l.upp.tri+m+m)])

    # SE parameters
    if (hess){
    SE.A<-diag(c(w$se[1:m]))
    SE.Chol<-diag(w$se[(m+1):(m*2)])
    nr.off.diag<-upper.tri(SE.Chol)
    l.upp.tri<-length(nr.off.diag[nr.off.diag==TRUE])
    
    SE.Chol[upper.tri(SE.Chol)] <- w$se[((length(diag(SE.A))*2)+1):((length(diag(SE.A))*2)+l.upp.tri)]
    SE.R<-t(SE.Chol)%*%SE.Chol
    
    SE.optima<-c(w$se[((length(diag(SE.A))*2)+l.upp.tri+1):((length(diag(SE.A))*2)+l.upp.tri+m)])
    
    SE.anc<-c(w$se[((length(diag(SE.A))*2)+l.upp.tri+1+m):((length(diag(SE.A))*2)+l.upp.tri+m+m)])
    }
    
    modelName<-"Multivariate model: OU; diagonal A matrix, symmetric R matrix"
    
  }

  if(A.matrix=="full" & R.matrix=="diag"){

    ### Eigenvalue decomposition of A ###
    A<-diag(c(w$par[1:m]))
    nr.off.diag<-upper.tri(A)
    l.upp.tri<-length(nr.off.diag[nr.off.diag==TRUE])
    A[upper.tri(A)] <- w$par[(length(diag(A))+1):(length(diag(A))+l.upp.tri)]
    A[lower.tri(A)] <- w$par[(length(diag(A))+l.upp.tri+1):(length(diag(A))+l.upp.tri+l.upp.tri)]


    ### The R (drift) matrix ###
    Chol<-diag(w$par[(length(diag(A))+l.upp.tri+l.upp.tri+1):(length(diag(A))+l.upp.tri+l.upp.tri+(length(diag(A))))])
    R<-t(Chol)%*%Chol

    ### Theta (optimal trait values) ###
    optima<-c(w$par[(length(diag(A))+l.upp.tri+l.upp.tri+(length(diag(A)))+1):(length(diag(A))+l.upp.tri+l.upp.tri+(length(diag(A)))+length(diag(A)))])

    ### The ancestral trait values ###
    ancestral.values<-c(w$par[(length(diag(A))+l.upp.tri+l.upp.tri+(length(diag(A)))+length(diag(A))+1):(length(diag(A))+l.upp.tri+l.upp.tri+(length(diag(A)))+length(diag(A))++length(diag(A)))])

    if (hess){
    # SE parameters
    SE.A<-diag(c(w$se[1:m]))
    nr.off.diag<-upper.tri(SE.A)
    l.upp.tri<-length(nr.off.diag[nr.off.diag==TRUE])
    SE.A[upper.tri(SE.A)] <- w$se[(length(diag(SE.A))+1):(length(diag(SE.A))+l.upp.tri)]
    SE.A[lower.tri(SE.A)] <- w$se[(length(diag(SE.A))+l.upp.tri+1):(length(diag(SE.A))+l.upp.tri+l.upp.tri)]
    
    
    ### The R (drift) matrix ###
    SE.Chol<-diag(w$se[(length(diag(SE.A))+l.upp.tri+l.upp.tri+1):(length(diag(SE.A))+l.upp.tri+l.upp.tri+(length(diag(SE.A))))])
    SE.R<-t(SE.Chol)%*%SE.Chol
    
    ### Theta (optimal trait values) ###
    SE.optima<-c(w$se[(length(diag(SE.A))+l.upp.tri+l.upp.tri+(length(diag(SE.A)))+1):(length(diag(SE.A))+l.upp.tri+l.upp.tri+(length(diag(SE.A)))+length(diag(SE.A)))])
    
    ### The ancestral trait values ###
    SE.anc<-c(w$se[(length(diag(SE.A))+l.upp.tri+l.upp.tri+(length(diag(SE.A)))+length(diag(SE.A))+1):(length(diag(SE.A))+l.upp.tri+l.upp.tri+(length(diag(SE.A)))+length(diag(SE.A))++length(diag(SE.A)))])
    }  
    
    modelName<-"Multivariate model: OU; full A matrix, diagonal R matrix"
    
  }

  if(A.matrix=="full" & R.matrix=="symmetric"){

    ### Eigenvalue decomposition of A ###
    A<-diag(c(w$par[1:m]))
    nr.off.diag<-upper.tri(A)
    l.upp.tri<-length(nr.off.diag[nr.off.diag==TRUE])
    A[upper.tri(A)] <- w$par[(length(diag(A))+1):(length(diag(A))+l.upp.tri)]
    A[lower.tri(A)] <- w$par[(length(diag(A))+l.upp.tri+1):(length(diag(A))+l.upp.tri+l.upp.tri)]

    ### The R (drift) matrix ###
    Chol<-diag(w$par[(length(diag(A))+l.upp.tri+l.upp.tri+1):(length(diag(A))+l.upp.tri+l.upp.tri+length(diag(A)))])
    nr.off.diag.R<-upper.tri(Chol)
    l.upp.tri.R<-length(nr.off.diag.R[nr.off.diag.R==TRUE])
    Chol[upper.tri(Chol)] <- w$par[(length(diag(A))+l.upp.tri+l.upp.tri+length(diag(A))+1):(length(diag(A))+l.upp.tri+l.upp.tri+length(diag(A))+l.upp.tri.R)]
    R<-t(Chol)%*%Chol

    ### Theta (optimal trait values) ###
    optima<-c(w$par[(length(diag(A))+l.upp.tri+l.upp.tri+length(diag(A))+l.upp.tri.R+1):(length(diag(A))+l.upp.tri+l.upp.tri+length(diag(A))+l.upp.tri.R+length(diag(A)))])

    ### The ancestral trait values ###
    ancestral.values<-c(w$par[(length(diag(A))+l.upp.tri+l.upp.tri+length(diag(A))+l.upp.tri.R+length(diag(A))+1):(length(diag(A))+l.upp.tri+l.upp.tri+length(diag(A))+l.upp.tri.R+length(diag(A))+length(diag(A)))])

    if (hess){
    # SE parameters
    SE.A<-diag(c(w$se[1:m]))
    nr.off.diag<-upper.tri(SE.A)
    l.upp.tri<-length(nr.off.diag[nr.off.diag==TRUE])
    SE.A[upper.tri(SE.A)] <- w$se[(length(diag(SE.A))+1):(length(diag(SE.A))+l.upp.tri)]
    SE.A[lower.tri(SE.A)] <- w$se[(length(diag(SE.A))+l.upp.tri+1):(length(diag(SE.A))+l.upp.tri+l.upp.tri)]
    
    ### The R (drift) matrix ###
    SE.Chol<-diag(w$se[(length(diag(SE.A))+l.upp.tri+l.upp.tri+1):(length(diag(SE.A))+l.upp.tri+l.upp.tri+length(diag(SE.A)))])
    nr.off.diag.R<-upper.tri(SE.Chol)
    l.upp.tri.R<-length(nr.off.diag.R[nr.off.diag.R==TRUE])
    SE.Chol[upper.tri(SE.Chol)] <- w$se[(length(diag(SE.A))+l.upp.tri+l.upp.tri+length(diag(SE.A))+1):(length(diag(SE.A))+l.upp.tri+l.upp.tri+length(diag(SE.A))+l.upp.tri.R)]
    SE.R<-t(SE.Chol)%*%SE.Chol
    
    ### Theta (optimal trait values) ###
    SE.optima<-c(w$se[(length(diag(SE.A))+l.upp.tri+l.upp.tri+length(diag(SE.A))+l.upp.tri.R+1):(length(diag(SE.A))+l.upp.tri+l.upp.tri+length(diag(SE.A))+l.upp.tri.R+length(diag(SE.A)))])
    
    ### The ancestral trait values ###
    SE.anc<-c(w$se[(length(diag(SE.A))+l.upp.tri+l.upp.tri+length(diag(SE.A))+l.upp.tri.R+length(diag(SE.A))+1):(length(diag(SE.A))+l.upp.tri+l.upp.tri+length(diag(SE.A))+l.upp.tri.R+length(diag(SE.A))+length(diag(SE.A)))])
    }
    
    modelName<-"Multivariate model: OU; full A matrix, symmetric R matrix"
    
  }


  if(A.matrix=="upper.tri" & R.matrix=="diag"){

    A<-diag(c(w$par[1:m]))
    nr.off.diag<-upper.tri(A)
    l.upp.tri<-length(nr.off.diag[nr.off.diag==TRUE])
    A[upper.tri(A)] <- w$par[(length(diag(A))+1):(length(diag(A))+l.upp.tri)]

    Chol<-diag(w$par[(length(diag(A))+l.upp.tri+1):((length(diag(A))+l.upp.tri+m))])
    R<-t(Chol)%*%Chol

    optima<-c(w$par[(length(diag(A))+l.upp.tri+m+1):((length(diag(A))+l.upp.tri+m+m))])

    ancestral.values<-c(w$par[(length(diag(A))+l.upp.tri+m+m+1):(length(diag(A))+l.upp.tri+m+m+m)])

    if (hess){
    # SE parameters
    SE.A<-diag(c(w$se[1:m]))
    nr.off.diag<-upper.tri(SE.A)
    l.upp.tri<-length(nr.off.diag[nr.off.diag==TRUE])
    SE.A[upper.tri(SE.A)] <- w$se[(length(diag(SE.A))+1):(length(diag(SE.A))+l.upp.tri)]
    
    SE.Chol<-diag(w$se[(length(diag(SE.A))+l.upp.tri+1):((length(diag(SE.A))+l.upp.tri+m))])
    SE.R<-t(SE.Chol)%*%SE.Chol
    
    SE.optima<-c(w$se[(length(diag(SE.A))+l.upp.tri+m+1):((length(diag(SE.A))+l.upp.tri+m+m))])
    
    SE.anc<-c(w$se[(length(diag(SE.A))+l.upp.tri+m+m+1):(length(diag(SE.A))+l.upp.tri+m+m+m)])
    }
    
    modelName<-"Multivariate model: OU; upper triangular A matrix, diagonal R matrix "
    
      }

  if(A.matrix=="upper.tri" & R.matrix=="symmetric"){

    A<-diag(c(w$par[1:m]))
    nr.off.diag.A<-upper.tri(A)
    l.upp.tri.A<-length(nr.off.diag.A[nr.off.diag.A==TRUE])
    A[upper.tri(A)] <- w$par[(length(diag(A))+1):(length(diag(A))+l.upp.tri.A)]

    Chol<-diag(w$par[(length(diag(A))+l.upp.tri.A+1):(length(diag(A))+l.upp.tri.A+m)])
    nr.off.diag.R<-upper.tri(Chol)
    l.upp.tri.R<-length(nr.off.diag.R[nr.off.diag.R==TRUE])
    Chol[upper.tri(Chol)] <- w$par[(length(diag(A))+l.upp.tri.A+m+1):(length(diag(A))+l.upp.tri.A+m+l.upp.tri.R)]
    R<-t(Chol)%*%Chol

    optima<-c(w$par[(length(diag(A))+l.upp.tri.A+m+l.upp.tri.R+1):(length(diag(A))+l.upp.tri.A+m+l.upp.tri.R+m)])

    ancestral.values<-c(w$par[(length(diag(A))+l.upp.tri.A+m+l.upp.tri.R+m+1):(length(diag(A))+l.upp.tri.A+m+l.upp.tri.R+m+m)])
    
    if (hess){
    # SE parameters
    SE.A<-diag(c(w$se[1:m]))
    nr.off.diag.A<-upper.tri(SE.A)
    l.upp.tri.A<-length(nr.off.diag.A[nr.off.diag.A==TRUE])
    SE.A[upper.tri(SE.A)] <- w$se[(length(diag(SE.A))+1):(length(diag(SE.A))+l.upp.tri.A)]
    
    SE.Chol<-diag(w$se[(length(diag(SE.A))+l.upp.tri.A+1):(length(diag(SE.A))+l.upp.tri.A+m)])
    nr.off.diag.R<-upper.tri(SE.Chol)
    l.upp.tri.R<-length(nr.off.diag.R[nr.off.diag.R==TRUE])
    SE.Chol[upper.tri(SE.Chol)] <- w$se[(length(diag(SE.A))+l.upp.tri.A+m+1):(length(diag(SE.A))+l.upp.tri.A+m+l.upp.tri.R)]
    SE.R<-t(SE.Chol)%*%SE.Chol
    
    SE.optima<-c(w$se[(length(diag(SE.A))+l.upp.tri.A+m+l.upp.tri.R+1):(length(diag(SE.A))+l.upp.tri.A+m+l.upp.tri.R+m)])
    
    SE.anc<-c(w$se[(length(diag(SE.A))+l.upp.tri.A+m+l.upp.tri.R+m+1):(length(diag(SE.A))+l.upp.tri.A+m+l.upp.tri.R+m+m)])
    }  
    
    modelName<-"Multivariate model: OU; upper triangular A matrix, symmetric R matrix "
    
  }

  if(A.matrix=="lower.tri" & R.matrix=="diag"){

    A<-diag(c(w$par[1:m]))
    nr.off.diag<-lower.tri(A)
    l.low.tri<-length(nr.off.diag[nr.off.diag==TRUE])
    A[lower.tri(A)] <- w$par[(length(diag(A))+1):(length(diag(A))+l.low.tri)]

    Chol<-diag(w$par[(length(diag(A))+l.low.tri+1):((length(diag(A))+l.low.tri+m))])
    R<-t(Chol)%*%Chol

    optima<-c(w$par[(length(diag(A))+l.low.tri+m+1):((length(diag(A))+l.low.tri+m+m))])

    ancestral.values<-c(w$par[(length(diag(A))+l.low.tri+m+m+1):(length(diag(A))+l.low.tri+m+m+m)])

    if (hess){
    # SE parameters
    SE.A<-diag(c(w$se[1:m]))
    nr.off.diag<-lower.tri(SE.A)
    l.low.tri<-length(nr.off.diag[nr.off.diag==TRUE])
    SE.A[lower.tri(SE.A)] <- w$se[(length(diag(SE.A))+1):(length(diag(SE.A))+l.low.tri)]
    
    SE.Chol<-diag(w$se[(length(diag(SE.A))+l.low.tri+1):((length(diag(SE.A))+l.low.tri+m))])
    SE.R<-t(SE.Chol)%*%SE.Chol
    
    SE.optima<-c(w$se[(length(diag(SE.A))+l.low.tri+m+1):((length(diag(SE.A))+l.low.tri+m+m))])
    
    SE.anc<-c(w$se[(length(diag(SE.A))+l.low.tri+m+m+1):(length(diag(SE.A))+l.low.tri+m+m+m)])
  }
  
    modelName<-"Multivariate model: OU; lower triangular A matrix, diagonal R matrix "
    
    }
  
  if(A.matrix=="lower.tri" & R.matrix=="symmetric"){

    A<-diag(c(w$par[1:m]))
    nr.off.diag.A<-lower.tri(A)
    l.low.tri.A<-length(nr.off.diag.A[nr.off.diag.A==TRUE])
    A[lower.tri(A)] <- w$par[(length(diag(A))+1):(length(diag(A))+l.low.tri.A)]

    Chol<-diag(w$par[(length(diag(A))+l.low.tri.A+1):(length(diag(A))+l.low.tri.A+m)])
    nr.off.diag.R<-upper.tri(Chol)
    l.upp.tri.R<-length(nr.off.diag.R[nr.off.diag.R==TRUE])
    Chol[upper.tri(Chol)] <- w$par[(length(diag(A))+l.low.tri.A+m+1):(length(diag(A))+l.low.tri.A+m+l.upp.tri.R)]
    R<-t(Chol)%*%Chol

    optima<-c(w$par[(length(diag(A))+l.low.tri.A+m+l.upp.tri.R+1):(length(diag(A))+l.low.tri.A+m+l.upp.tri.R+m)])

    ancestral.values<-c(w$par[(length(diag(A))+l.low.tri.A+m+l.upp.tri.R+m+1):(length(diag(A))+l.low.tri.A+m+l.upp.tri.R+m+m)])
 
    if (hess){
    # SE parameters
    SE.A<-diag(c(w$se[1:m]))
    nr.off.diag.A<-lower.tri(SE.A)
    l.low.tri.A<-length(nr.off.diag.A[nr.off.diag.A==TRUE])
    SE.A[lower.tri(SE.A)] <- w$se[(length(diag(SE.A))+1):(length(diag(SE.A))+l.low.tri.A)]
    
    SE.Chol<-diag(w$se[(length(diag(SE.A))+l.low.tri.A+1):(length(diag(SE.A))+l.low.tri.A+m)])
    nr.off.diag.R<-upper.tri(SE.Chol)
    l.upp.tri.R<-length(nr.off.diag.R[nr.off.diag.R==TRUE])
    SE.Chol[upper.tri(SE.Chol)] <- w$se[(length(diag(SE.A))+l.low.tri.A+m+1):(length(diag(SE.A))+l.low.tri.A+m+l.upp.tri.R)]
    SE.R<-t(SE.Chol)%*%SE.Chol
    
    SE.optima<-c(w$se[(length(diag(SE.A))+l.low.tri.A+m+l.upp.tri.R+1):(length(diag(SE.A))+l.low.tri.A+m+l.upp.tri.R+m)])
    
    SE.anc<-c(w$se[(length(diag(SE.A))+l.low.tri.A+m+l.upp.tri.R+m+1):(length(diag(SE.A))+l.low.tri.A+m+l.upp.tri.R+m+m)])
    }
    
    modelName<-"Multivariate model: OU; lower triangular A matrix, symmetric R matrix "
    
     }

  if(A.matrix=="OUBM"){

    A<-diag(c(w$par[1:(m-1)],0))
    nr.off.diag.A<-upper.tri(A)
    l.upp.tri.A<-length(nr.off.diag.A[nr.off.diag.A==TRUE])
    A[upper.tri(A)] <- w$par[(m):(m-1+l.upp.tri.A)]

    Chol<-diag(w$par[(m+l.upp.tri.A):(m+m+l.upp.tri.A-1)])
    R<-t(Chol)%*%Chol

    optima<-c(w$par[(m+m+l.upp.tri.A):(m+m+m+l.upp.tri.A-1)])
    optima[length(optima)]<-NA

    ancestral.values<-c(w$par[(m+m+m+l.upp.tri.A):(m+m+m+m+l.upp.tri.A-1)])
    K <- length(init.par)-1
    
    if (hess){
    # SE parameters
    SE.A<-diag(c(w$se[1:(m-1)],0))
    nr.off.diag.A<-upper.tri(SE.A)
    l.upp.tri.A<-length(nr.off.diag.A[nr.off.diag.A==TRUE])
    SE.A[upper.tri(SE.A)] <- w$se[(m):(m-1+l.upp.tri.A)]
    
    SE.Chol<-diag(w$se[(m+l.upp.tri.A):(m+m+l.upp.tri.A-1)])
    SE.R<-t(SE.Chol)%*%SE.Chol
    
    SE.optima<-c(w$se[(m+m+l.upp.tri.A):(m+m+m+l.upp.tri.A-1)])
    SE.optima[length(SE.optima)]<-NA
    
    SE.anc<-c(w$se[(m+m+m+l.upp.tri.A):(m+m+m+m+l.upp.tri.A-1)])
    }
    
    modelName<-"Multivariate model: OUBM"
      }

  half.life<-log(2)/diag(A)
  
  

  if (hess){
    wc<-as.evoTS.multi.OU.fit(converge, modelName = modelName, logL = w$value, ancestral.values = ancestral.values, SE.anc = SE.anc, optima = optima, SE.optima = SE.optima, A = A, SE.A = SE.A, half.life = half.life, R = R, SE.R = SE.R,
                                                       method = "Joint", K = K, n = length(yy$xx[,1]), iter=iter)
  }
  
  if (hess == FALSE){
    SE.anc <- NA
    SE.optima <- NA
    SE.A <- NA 
    SE.R <- NA 
    wc<-as.evoTS.multi.OU.fit(converge, modelName = modelName, logL = w$value, ancestral.values = ancestral.values, SE.anc = SE.anc, optima = optima, SE.optima = SE.optima, A = A, SE.A = SE.A, half.life = half.life, R = R, SE.R = SE.R,
                                 method = "Joint", K = K, n = length(yy$xx[,1]), iter=iter)
  }
  return(wc)

 }
