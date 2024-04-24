#' @title Fit separate multivariate Unbiased Random Walk models to two different segments of a multivariate evolutionary sequence (time-series).
#'
#' @description Function to find maximum likelihood solutions for multivariate Unbiased Random Walk models fitted to two different segments of a multivariate evolutionary sequence (time-series).
#'
#' @param yy a multivariate evoTS object.
#'
#' @param minb the minimum number of samples within a segment to consider.
#'
#' @param hess logical, indicating whether to calculate standard errors from the Hessian matrix.
#' 
#' @param pool indicating whether to pool variances across samples
#'
#' @param shift.point the sample in the time series that represents the first sample in the second segment.
#'
#' @param method optimization method, passed to function optim. Default is "L-BFGS-B".
#'
#' @param trace logical, indicating whether information on the progress of the optimization is printed.
#'
#' @param iterations the number of times the optimization method is run from different starting points. Default is NULL, meaning the optimization is run once.
#'
#' @param iter.sd defines the standard deviation of the Gaussian distribution from which starting values for the optimization routine is run. Default is 1.
#'
#' @details The function searches - using an optimization routine - for the maximum-likelihood solution for a multivariate Unbiased Random Walk model ti two non-overlapping segments in the time series.
#'
#' The argument 'method' is passed to the 'optim' function and is included for the convenience of users to better control the optimization routine. The the default method (L-BFGS-B) seems to work for most evolutionary sequences.
#'
#' Initial estimates to start the optimization come from maximum-likelihood estimates of the univariate Unbiased Random Walk model (from the paleoTS package) fitted to each time-series separately.
#'
#' It is good practice to repeat any numerical optimization procedure from different starting points. This is especially important for complex models as the log-likelihood surface might contain more than one peak. The number of iterations is controlled by the argument 'iterations'. The function will report the model parameters from the iteration with the highest log-likelihood.
#'
#'@return First part of the output reports the log-likelihood of the model and its AICc score. The second part of the output is the maximum log-likelihood model parameters (ancestral.values, R). The last part of the output gives information about the number of parameters in the model (K), number of samples in the data (n) and number of times the optimization routine was run (iter).
#'
#'@note The models have been implemented to be compatible with the joint parameterization routine in the package paleoTS. The optimization is therefore fit using the actual sample values, with the autocorrelation among samples accounted for in the log-likelihood function. The joint distribution of sample means is multivariate normal, with means and variance-covariances determined by evolutionary parameters and sampling errors.
#'
#'@author Kjetil Lysne Voje
#'
#'@references Revell, L. J. & Harmon, L. Testing quantitative genetic hypotheses about the evolutionary rate matrix for continuous characters. \emph{Evolutionary Ecology Research} 10, 311–331 (2008).
#'
#'@export
#'
#'@examples
#'## Generate an evoTS object by simulating a multivariate dataset
#'x <- sim.multi.URW(60)
#'
#'## Fit two multivariate Unbiased Random Walk models to separate parts of the time-series.
#'fit.multivariate.URW.shift(x, shift.point = 31)

fit.multivariate.URW.shift<-function (yy, minb = 10, hess = FALSE, pool = TRUE, shift.point = NULL, method = "L-BFGS-B", trace=FALSE, iterations=NULL, iter.sd=NULL)
{
  
  ng<-2 # The program is currently constrained to fitting two R matrices.
  ns <- nrow(yy$xx)
  m <- ncol(yy$xx) # number of traits
  
  if (pool==TRUE) { 
    for (i in 1:m){
      
      tmp<-paleoTS::as.paleoTS(yy$xx[,i], yy$vv[,i], yy$nn[,i], yy$tt[,i])
      tmp<- paleoTS::pool.var(tmp, ret.paleoTS = TRUE)
      yy$vv[,i]<-tmp$vv
    }
  }
  
  {
    if (is.null(shift.point) == FALSE & length(shift.point) > 1) {stop("A user-defined shift point can only be a single shift. Make sure to only define a single shift point when using the shift.point argument")}
  }
  
   #### If a switch point has been defined and iterations = 0 (NULL)
 
   ##### Start of non-iteration routine #####
  if(is.null(shift.point) == FALSE) {
  
  GG <-shift.point
  GG<-as.matrix(GG)
  
  #Define number of shift points:
  nc <- ncol(GG)
  
  n <- nrow(yy$xx) # number of samples/populations
  gg <- shift2gg(GG[, 1], ns)
  states <- length(unique(gg))

  #define vectors based on length of first and second segments in the time series
  seg.1 <- which(gg == 1)
  seg.2 <- which(gg == 2)

  #Extract the time vector for each random walks:
  tt.seg.1 <- yy$tt[,1][seg.1]
  tt.seg.2 <- yy$tt[,1][seg.2]

  #Create distance matrices
  C.1<-matrix(0, n, n)
  C.1[seg.1, seg.1]<- outer(tt.seg.1, tt.seg.1, FUN = pmin)
  C.2<-matrix(0, n, n)
  C.2 [seg.2, seg.2] <- outer(tt.seg.2, tt.seg.2, FUN = pmin)

  C <- list(C.1 = C.1,  C.2= C.2)

  X <- yy$xx # Character matrix with dimensions n * m
  y <- as.matrix(as.vector(X)) # Vectorized version of X


  ### Estimate two full R matrices

  # Define initial parameter values for the optimization routine
  init.trait.var<-apply(yy$xx,2,var)
  temp.matrix<-cov(as.matrix(yy$xx))
  init.cov.traits<-unique(temp.matrix[row(temp.matrix)!=col(temp.matrix)])
  anc.values<-yy$xx[1,]

  init.par<-c(init.trait.var, init.trait.var, init.cov.traits, init.cov.traits, anc.values)
  lower.limit<-c(rep(0,length(init.trait.var)), rep(0,length(init.trait.var)),  rep(NA,length(init.cov.traits)), rep(NA, length(init.cov.traits)), rep(NA, length(anc.values)))

  if (is.numeric(iterations) == FALSE){
  
   if (method == "L-BFGS-B")  {
  w<-optim(init.par, fn = logL.joint.multi.R, C = C, y = y, m = m, n = n, anc.values = anc.values,  yy = yy,
             control = list(fnscale = -1, maxit=10000, trace = trace), method = "L-BFGS-B", hessian = hess, lower = lower.limit)
                              }

  if (method == "Nelder-Mead")  {
  w<-optim(init.par, fn = logL.joint.multi.R, C = C, y = y, m = m, n = n, anc.values = anc.values,  yy = yy,
             control = list(fnscale = -1, maxit=10000, trace = trace), method = "Nelder-Mead" , hessian = hess)
  }
  if (method == "SANN")  {
    w<-optim(init.par, fn = logL.joint.multi.R, C = C, y = y, m = m, n = n, anc.values = anc.values, yy = yy,
             control = list(fnscale = -1, maxit=10000, trace = trace), method = "SANN" , hessian = hess, lower = lower.limit)
  }
    
  }

  else
  {

    #### If a switch point has been defined and iterations have been defined:
    
    ### Start iterations from different starting values
      if(is.numeric(iter.sd) == FALSE) iter.sd <-1
      #if(is.numeric(max.attemps) == FALSE) max.attemps <-100000
      log.lik.tmp<-rep(NA, 1000000)
      www<-list()
      
      for (k in 1:1000000){
        tryCatch({
          init.par_temp<-init.par
          init.par<-rnorm(length(init.par_temp), init.par_temp, iter.sd)
          
          
          if (method == "L-BFGS-B")  {
            www[[k]]<-optim(init.par, fn = logL.joint.multi.R, C = C, y = y, m = m, n = n, anc.values = anc.values,  yy = yy,
                            control = list(fnscale = -1, maxit=10000, trace = trace), method = "L-BFGS-B", hessian = hess, lower = lower.limit)
          }
          
          if (method == "Nelder-Mead")  {
            www[[k]]<-optim(init.par, fn = logL.joint.multi.R, C = C, y = y, m = m, n = n, anc.values = anc.values,  yy = yy,
                            control = list(fnscale = -1, maxit=10000, trace = trace), method = "Nelder-Mead" , hessian = hess)
          }
          if (method == "SANN")  {
            www[[k]]<-optim(init.par, fn = logL.joint.multi.R, C = C, y = y, m = m, n = n, anc.values = anc.values,  yy = yy,
                            control = list(fnscale = -1, maxit=10000, trace = trace), method = "SANN" , hessian = hess, lower = lower.limit)
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
        if(max(na.exclude(log.lik.tmp)) == www_reduced[[j]]$value) w<-www_reduced[[j]]
      }
  }
  
  }
  
  
  #### If a switch point has NOT been defined and iterations = 0 (NULL)
  
  ##### Start of non-iteration routine #####

  
  if (is.numeric(shift.point) == FALSE) { #Slutter ved 330
    {
      if (ns/2 < minb) {stop("Number of samples in time-series is too small compared to the minimum number of samples within each. Try setting minb to a smaller number")}
    }

    GG <- shifts(ns, ng, minb = minb)
    GG<-as.matrix(GG)
    
    #Define number of shift points:
    nc <- ncol(GG)
    
    cat("Total # hypotheses: ", nc, "\n")
    #Create empty list
    wl <- list()
   # best.run <- list()
    #create array with length = to shift points and where every entry = -Inf
    logl <- array(-Inf, dim = nc)
    
    #start loop for estimating maximum likelihood parameters for each data set defined by the shift points
    for (i in 1:nc) {
      cat(i, " ")
      n <- nrow(yy$xx) # number of samples/populations
      gg <- shift2gg(GG[, i], ns)
      states <- length(unique(gg))
      
      #define vectors based on length of first and second segments in the time series
      seg.1 <- which(gg == 1)
      seg.2 <- which(gg == 2)
      
      #Extract the time vector for each random walks:
      tt.seg.1 <- yy$tt[,1][seg.1]
      tt.seg.2 <- yy$tt[,1][seg.2]
      
      #Create distance matrices
      C.1<-matrix(0, n, n)
      C.1[seg.1, seg.1]<- outer(tt.seg.1, tt.seg.1, FUN = pmin)
      C.2<-matrix(0, n, n)
      C.2 [seg.2, seg.2] <- outer(tt.seg.2, tt.seg.2, FUN = pmin)
      
      C <- list(C.1 = C.1,  C.2= C.2)
      
      X <- yy$xx # Character matrix with dimensions n * m
      y <- as.matrix(as.vector(X)) # Vectorized version of X
      
      
      ### Estimate two full R matrices
      
      # Define initial parameter values for the optimization routine
      init.trait.var<-apply(yy$xx,2,var)
      temp.matrix<-cov(as.matrix(yy$xx))
      init.cov.traits<-unique(temp.matrix[row(temp.matrix)!=col(temp.matrix)])
      anc.values<-yy$xx[1,]
      
      init.par<-c(init.trait.var, init.trait.var, init.cov.traits, init.cov.traits, anc.values)
      lower.limit<-c(rep(0,length(init.trait.var)), rep(0,length(init.trait.var)),  rep(NA,length(init.cov.traits)), rep(NA, length(init.cov.traits)), rep(NA, length(anc.values)))
      
      if(is.numeric(iterations) == FALSE) { 
      
      if (method == "L-BFGS-B")  {
        w<-optim(init.par, fn = logL.joint.multi.R, C = C, y = y, m = m, n = n, anc.values = anc.values,  yy = yy,
                 control = list(fnscale = -1, maxit=10000, trace = trace), method = "L-BFGS-B", hessian = hess, lower = lower.limit)
      }
      
      if (method == "Nelder-Mead")  {
        w<-optim(init.par, fn = logL.joint.multi.R, C = C, y = y, m = m, n = n, anc.values = anc.values,  yy = yy,
                 control = list(fnscale = -1, maxit=10000, trace = trace), method = "Nelder-Mead" , hessian = hess)
      }
      if (method == "SANN")  {
        w<-optim(init.par, fn = logL.joint.multi.R, C = C, y = y, m = m, n = n, anc.values = anc.values, yy = yy,
                 control = list(fnscale = -1, maxit=10000, trace = trace), method = "SANN" , hessian = hess, lower = lower.limit)
      }
      
      logl[i] <- w$value
      wl[[i]] <- w
      }
      
      else {
        
        if(is.numeric(iter.sd) == FALSE) iter.sd <-1
        #if(is.numeric(max.attemps) == FALSE) max.attemps <-100000
        log.lik.tmp<-rep(NA, 1000000)
        www<-list()
        
        for (k in 1:1000000){
          tryCatch({
            init.par_temp<-init.par
            init.par<-rnorm(length(init.par_temp), init.par_temp, iter.sd)
            
            
            if (method == "L-BFGS-B")  {
              www[[k]]<-optim(init.par, fn = logL.joint.multi.R, C = C, y = y, m = m, n = n, anc.values = anc.values,  yy = yy,
                              control = list(fnscale = -1, maxit=10000, trace = trace), method = "L-BFGS-B", hessian = hess, lower = lower.limit)
            }
            
            if (method == "Nelder-Mead")  {
              www[[k]]<-optim(init.par, fn = logL.joint.multi.R, C = C, y = y, m = m, n = n, anc.values = anc.values,  yy = yy,
                              control = list(fnscale = -1, maxit=10000, trace = trace), method = "Nelder-Mead" , hessian = hess)
            }
            if (method == "SANN")  {
              www[[k]]<-optim(init.par, fn = logL.joint.multi.R, C = C, y = y, m = m, n = n, anc.values = anc.values,  yy = yy,
                              control = list(fnscale = -1, maxit=10000, trace = trace), method = "SANN" , hessian = hess, lower = lower.limit)
            }
            log.lik.tmp[k]<-www[[k]]$value
          }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
          
          if (length(na.exclude(log.lik.tmp)) == iterations){
            break
          }
        }
        
        # Need to remove entries in www in case there are iterations where the initial parameter estimates did not work.
        www_tmp<-list()
        for (o in 1:k){
          if (is.character(www[[o]][1]) == FALSE) www_tmp[o]<-list(www[[o]])
        }
        
        www_reduced<-www_tmp[!sapply(www_tmp,is.null)]
        
        
        for (j in 1:length(www_reduced)){
          if(max(na.exclude(log.lik.tmp)) == www_reduced[[j]]$value) wl[[i]]<-www_reduced[[j]]
        }
        
       logl[i] <- wl[[i]]$value
        #wl[[i]] <- best.run  
        
      } #else
      #logl[i]<-unlist(lapply(best.run, "[[", 2))
      
    }   #for each switch point 
    
  #  for (u in 1:nc){
  #    if(max(logl) == best.run[[u]]$value) wl<-best.run[[u]]
  #  }
    
    #wl[[i]] <- best.run
    #print(logl)
    
    winner <- which.max(logl)
    w <- wl[[winner]]
    ss <- GG[, winner]
    names(ss) <- paste("shift", 1:(ng - 1), sep = "")
    w$par <- append(w$par, ss)
    w$all.logl <- logl
    w$GG <- GG
  }
  
  
  
    # number of parameters
    K <- length(init.par) + (states-1) #parameters in the R matrices + ancestral values for each trait + number of shift points


  if (is.numeric(iterations) ==TRUE) {
  iter<-iterations
 # w<-best.run
  }
    
    if (w$convergence == 10) converge<-"The search algorithm stopped as it did not make progress towards the optimal solution"
    if (w$convergence == 0) converge<-"Model converged successfully"
    if (w$convergence == 1) converge<-"The maximum number of iterations was reached and the search algorithm exited"
    if (w$convergence == 51) converge<-"The model did not converge due to en error in L-BFGS-B. Reported estimates are not the maximum likelihood"
    if (w$convergence == 52) converge<-"The model did not converge due to en error in L-BFGS-B. Reported estimates are not the maximum likelihood"

  if (hess) {
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
    SE.R1<-matrix(0, nrow=m, ncol=m)
    SE.R2<-matrix(0, nrow=m, ncol=m)
    diag(SE.R1)<-w$se[1:m]
    diag(SE.R2)<-w$se[(m+1):(m+m)]

    locations.SE.R1<-which(SE.R1 == 0, arr.ind = T)
    location.upper.tri.R<-which(locations.SE.R1[,1] < locations.SE.R1[,2])

    upper.first<-w$se[(m+m+1):(m+m+length(location.upper.tri.R))]
    upper.second<-w$se[(m+m+length(location.upper.tri.R)+1):(m+m+length(location.upper.tri.R)+length(location.upper.tri.R))]

    for (i in 1:m){
      SE.R1[locations.SE.R1[,1][location.upper.tri.R[i]],locations.SE.R1[,2][location.upper.tri.R[i]]]<-upper.first[i]
      SE.R2[locations.SE.R1[,1][location.upper.tri.R[i]],locations.SE.R1[,2][location.upper.tri.R[i]]]<-upper.second[i]
    }

    SE.R.1<-t(SE.R1)%*%SE.R1
    SE.R.2<-t(SE.R2)%*%SE.R2
    SE.R <- list(SE.R.1 = SE.R.1,  SE.R.2 = SE.R.2)
    SE.anc <- w$se[((length(init.trait.var)*2 + length(init.cov.traits)*2)+1) : length(init.par)]
  }
  if (hess == FALSE) {
    SE.R<-NA
    SE.anc <-NA
  }

  chole.1<-matrix(0, nrow=m, ncol=m)
  chole.2<-matrix(0, nrow=m, ncol=m)
  diag(chole.1)<-w$par[1:m]
  diag(chole.2)<-w$par[(m+1):(m+m)]

  locations.R1<-which(chole.1 == 0, arr.ind = T)
  location.upper.tri.R<-which(locations.R1[,1] < locations.R1[,2])

  upper.first<-w$par[(m+m+1):(m+m+length(location.upper.tri.R))]
  upper.second<-w$par[(m+m+length(location.upper.tri.R)+1):(m+m+length(location.upper.tri.R)+length(location.upper.tri.R))]

  for (i in 1:m){
    chole.1[locations.R1[,1][location.upper.tri.R[i]],locations.R1[,2][location.upper.tri.R[i]]]<-upper.first[i]
    chole.2[locations.R1[,1][location.upper.tri.R[i]],locations.R1[,2][location.upper.tri.R[i]]]<-upper.second[i]
  }

  R.1<-t(chole.1)%*%chole.1
  R.2<-t(chole.2)%*%chole.2

  R <- list(R.1 = R.1,  R.2= R.2)

  if (is.numeric(iterations) == FALSE) {
    iter<-NA
  }

  ancestral.values<-w$par[((length(init.trait.var)*2 + length(init.cov.traits)*2)+1) : length(init.par)]
  
  if(is.null(shift.point) == TRUE) shift.point<-unname(w$par[length(init.par)+1]) 
  
  if(is.null(shift.point) == FALSE) shift.point<-unname(shift.point) 
  

  wc<-as.evoTS.multi.URW.shift.fit(converge, modelName = "Multivariate model: Random walk with two R matrices (with non-zero off-diagonal elements)", logL = w$value, ancestral.values = ancestral.values, SE.anc = SE.anc, R = R, SE.R = SE.R,
                            shift.point = shift.point, method = "Joint", K = K, n = length(yy$xx[,1]), iter=iter)
  return(wc)
}
