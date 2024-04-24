#' @title Fit two models to two separate segments to an evolutionary sequence (time-series).
#'
#' @description Wrapper function to find maximum likelihood solutions to two models to an evolutionary sequence.
#'
#' @param y an univariate evoTS object.
#'
#'@param model1 the model fitted to the first segment. Options are Stasis, URW, GRW, OU.
#'
#'@param model2 the model fitted to the second segment. Options are Stasis, URW, GRW, OU.
#'
#'@param fit.all logical indicating whether to fit all pairwise combinations of the four models to the evolutionary sequence (time-series).
#'
#'@param minb the minimum number of samples within a segment to consider
#'
#'@param shift.point The sample that split the time-series into two segments. The samples are passed to the argument as a vector. Default is NULL, which means all possible shift points will be assessed constrained by how minb is defined.
#'
#' @param pool logical indicating whether to pool variances across samples
#'
#' @param silent if TRUE, less information is printed to the screen as the model is fit
#'
#' @param hess logical, indicating whether to calculate standard errors from the Hessian matrix.
#'
#'#'
#'@return the function  returns a list of all investigated models and their highest log-likelihood (and their corresponding AICc and AICc weight).
#'\item{logL}{the log-likelihood of the optimal solution}
#'\item{AICc}{AIC with a correction for small sample sizes}
#'\item{parameters}{parameter estimates}
#'\item{modelName}{abbreviated model name}
#'\item{method}{Joint consideration of all samples}
#'\item{K}{number of parameters in the model}
#'\item{n}{the number of observations/samples}
#'\item{all.logl}{log-likelihoods for all tested partitions of the series into segments. Will return a single value if shift points have been given}
#'\item{GG}{matrix of indices of initial samples of each tested segment configuration; each column of GG corresponds to the elements of all.logl}
#'
#'In addition, if fit.all=TRUE the function also returns a list of all investigated models and their highest log-likelihood (and their corresponding AICc and AICc weight).
#'
#'@note The models have been implemented to be compatible with the joint parameterization routine in the package paleoTS. The optimization is therefore fit using the actual sample values, with the autocorrelation among samples accounted for in the log-likelihood function. The joint distribution of sample means is multivariate normal, with means and variance-covariances determined by evolutionary parameters and sampling errors.
#'
#'@author Kjetil Lysne Voje
#'
#'@references Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond. \emph{Paleobiology} 32:578–601
#'@references Hunt, G., Bell, M. A. & Travis, M. P. Evolution towards a new adaptive optimum: Phenotypic evolution in a fossil stickleback lineage. \emph{Evolution} 62:700–710 (2008)
#'
#'@export
#'
#'@examples
#'
#'##Generate a paleoTS object.
#'x <- paleoTS::sim.GRW(30)
#'
#'## Fit a mode-shift model without defining a shift point (the example may take > 5 seconds to run)
#'fit.mode.shift(x, model1="URW", model2="Stasis")


fit.mode.shift<-function (y, model1=c("Stasis", "URW", "GRW", "OU"), model2=c("Stasis", "URW", "GRW", "OU"), fit.all=FALSE, minb = 7, shift.point = NULL, pool = TRUE, silent = FALSE, hess = FALSE)
{
  ns <- length(y$mm)
  ng <- 2
  y$start.age<-NULL
  if(is.numeric(shift.point) == TRUE) GG <-shift.point else GG <- shifts(ns, ng, minb = minb)
  GG<-as.matrix(GG)
  if (ncol(GG) == 1) print("Fitting the model for a user-defined shift point") else print("Searching all possible shift points in the evolutionary sequence")


  if(fit.all == TRUE){
  #Define number of shift points:
  nc <- ncol(GG)
  #Create empty list
  wl <- list()
  #create array with length = to shift points and where every entry ) -Inf
  logl <- array(-Inf, dim = nc)

  a1 <- vector(mode = "list", length = nc)
  a2 <- vector(mode = "list", length = nc)
  a3 <- vector(mode = "list", length = nc)
  a4 <- vector(mode = "list", length = nc)
  a5 <- vector(mode = "list", length = nc)
  a6 <- vector(mode = "list", length = nc)
  a7 <- vector(mode = "list", length = nc)
  a8 <- vector(mode = "list", length = nc)
  a9 <- vector(mode = "list", length = nc)
  a10 <- vector(mode = "list", length = nc)
  a11 <- vector(mode = "list", length = nc)
  a12 <- vector(mode = "list", length = nc)
  a13 <- vector(mode = "list", length = nc)
  a14 <- vector(mode = "list", length = nc)
  a15 <- vector(mode = "list", length = nc)
  a16 <- vector(mode = "list", length = nc)

  for (i in 1:nc) {
    if (!silent)
      cat(i, " ")
    #defines which data point in the time series that belong to each of the two sets
    gg <- shift2gg(GG[, i], ns)

  ## Her må a1 være en liste med lengde ns, slik at alle resultatene kan lagres for alle i (shift points)
    a1[[i]]  <- opt.joint.Stasis.Stasis(y, gg, pool = pool, hess = hess, meth = "L-BFGS-B") #Done
    a2[[i]]  <- opt.joint.Stasis.URW(y, gg, pool = pool, hess = hess, meth = "L-BFGS-B") #Done
    a3[[i]]  <- opt.joint.Stasis.GRW(y, gg, pool = pool, hess = hess, meth = "L-BFGS-B") #Done
    a4[[i]]  <- opt.joint.Stasis.OU(y, gg, pool = pool, hess = hess, meth = "L-BFGS-B") #Done
    a5[[i]]  <- opt.joint.URW.URW(y, gg, pool = pool, hess = hess, meth = "L-BFGS-B") #Done
    a6[[i]]  <- opt.joint.URW.GRW(y, gg, pool = pool, hess = hess, meth = "L-BFGS-B") #Done
    a7[[i]]  <- opt.joint.URW.OU(y, gg, pool = pool, hess = hess, meth = "L-BFGS-B") #Done
    a8[[i]]  <- opt.joint.GRW.GRW(y, gg, pool = pool, hess = hess, meth = "L-BFGS-B") #Done
    a9[[i]]  <- opt.joint.GRW.OU(y, gg, pool = pool, hess = hess, meth = "L-BFGS-B") # Done
    a10[[i]] <- opt.joint.OU.OU(y, gg, pool = pool, hess = hess, meth = "L-BFGS-B") #Done
    a11[[i]] <- opt.joint.OU.GRW(y, gg, pool = pool, hess = hess, meth = "L-BFGS-B") #Done
    a12[[i]] <- opt.joint.OU.URW(y, gg, pool = pool, hess = hess, meth = "L-BFGS-B") #Done
    a13[[i]] <- opt.joint.OU.Stasis(y, gg, pool = pool, hess = hess, meth = "L-BFGS-B") #Done
    a14[[i]] <- opt.joint.GRW.URW(y, gg, pool = pool, hess = hess, meth = "L-BFGS-B")
    a15[[i]] <- opt.joint.GRW.Stasis(y, gg, pool = pool, hess = hess, meth = "L-BFGS-B")
    a16[[i]] <- opt.joint.URW.Stasis(y, gg, pool = pool, hess = hess, meth = "L-BFGS-B")
  }
    wl<-list(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16)
    # Plukke ut hvilken av de 16 modellene som er best, for så å sammenlikne dem.
winner<-rep(NA, 16)
for (j in 1:16){
  tmp<-rep(NA, nc)
    for (i in 1:nc){
      tmp[i]<-wl[[j]][[i]]$AICc
    }
winner[j]<-which(tmp == min(tmp))
}

best.model <- which.min(c(wl[[1]][[winner[1]]]$AICc, wl[[2]][[winner[2]]]$AICc,wl[[3]][[winner[3]]]$AICc, wl[[4]][[winner[4]]]$AICc, wl[[5]][[winner[5]]]$AICc,
                   wl[[6]][[winner[6]]]$AICc, wl[[7]][[winner[7]]]$AICc,wl[[8]][[winner[8]]]$AICc, wl[[9]][[winner[9]]]$AICc, wl[[10]][[winner[10]]]$AICc,
                   wl[[11]][[winner[11]]]$AICc, wl[[12]][[winner[12]]]$AICc,wl[[13]][[winner[13]]]$AICc, wl[[14]][[winner[14]]]$AICc, wl[[15]][[winner[15]]]$AICc,
                   wl[[16]][[winner[16]]]$AICc))

#winner <- wl[tmp]
ww <- wl[[best.model]][[winner[best.model]]]
iteration<-winner[best.model]
ss <- GG[, iteration]
names(ss) <- paste("shift", 1:(ng - 1), sep = "")
ww$parameters <- append(ww$parameters, ss)
ww$GG <- ss
mc<-paleoTS::compareModels(wl[[1]][[winner[1]]], wl[[2]][[winner[2]]],wl[[3]][[winner[3]]], wl[[4]][[winner[4]]], wl[[5]][[winner[5]]],
              wl[[6]][[winner[6]]], wl[[7]][[winner[7]]],wl[[8]][[winner[8]]], wl[[9]][[winner[9]]], wl[[10]][[winner[10]]],
              wl[[11]][[winner[11]]], wl[[12]][[winner[12]]],wl[[13]][[winner[13]]], wl[[14]][[winner[14]]], wl[[15]][[winner[15]]],
              wl[[16]][[winner[16]]], silent = FALSE)

invisible(mc)
    out<-list(ww)

  }
  else {

      #Define number of shift points:
      nc <- ncol(GG)
      if (!silent)
        cat("Total # hypotheses: ", nc, "\n")
      #Create empty list
      wl <- list()
      #create array with length = to shift points and where every entry ) -Inf
      logl <- array(-Inf, dim = nc)

      for (i in 1:nc) {
        if (!silent)
          cat(i, " ")
        #defines which data point in the time series that belong to each of the two sets
        gg <- shift2gg(GG[, i], ns)

  if (model1 == "Stasis" & model2=="Stasis"){
    w <- opt.joint.Stasis.Stasis(y, gg, pool = pool, hess = hess) #Done
  }
    if (model1 == "Stasis" & model2=="URW"){
      w <- opt.joint.Stasis.URW(y, gg, pool = pool, hess = hess) #Done
    }
    if (model1 == "Stasis" & model2=="GRW"){
      w <- opt.joint.Stasis.GRW(y, gg, pool = pool, hess = hess) #Done
    }
    if (model1 == "Stasis" & model2=="OU"){
      w <- opt.joint.Stasis.OU(y, gg, pool = pool, hess = hess) #Done
    }
    if (model1 == "URW" & model2=="URW"){
      w <- opt.joint.URW.URW(y, gg, pool = pool, hess = hess) #Done
    }
    if (model1 == "URW" & model2=="GRW"){
      w <- opt.joint.URW.GRW(y, gg, pool = pool, hess = hess) #Done
    }
    if (model1 == "URW" & model2=="OU"){
      w <- opt.joint.URW.OU(y, gg, pool = pool, hess = hess) #Done
    }
    if (model1 == "GRW" & model2=="GRW"){
      w <- opt.joint.GRW.GRW(y, gg, pool = pool, hess = hess) #Done
    }
    if (model1 == "GRW" & model2=="OU"){
      w <- opt.joint.GRW.OU(y, gg, pool = pool, hess = hess) # Done
    }
    if (model1 == "OU" & model2=="OU"){
      w <- opt.joint.OU.OU(y, gg, pool = pool, hess = hess) #Done
    }
    if (model1 == "OU" & model2=="GRW"){
      w <- opt.joint.OU.GRW(y, gg, pool = pool, hess = hess) #Done
    }
    if (model1 == "OU" & model2=="URW"){
      w <- opt.joint.OU.URW(y, gg, pool = pool, hess = hess) #Done
    }
    if (model1 == "OU" & model2=="Stasis"){
      w <- opt.joint.OU.Stasis(y, gg, pool = pool, hess = hess) #Done
    }
    if (model1 == "GRW" & model2=="URW"){ #Done
      w <- opt.joint.GRW.URW(y, gg, pool = pool, hess = hess)
    }
    if (model1 == "GRW" & model2=="Stasis"){ #Done
      w <- opt.joint.GRW.Stasis(y, gg, pool = pool, hess = hess)
    }
    if (model1 == "URW" & model2=="Stasis"){
      w <- opt.joint.URW.Stasis(y, gg, pool = pool, hess = hess)
    }
    logl[i] <- w$logL
    wl[[i]] <- w
  }

  if (!silent)
    cat("\n")
  winner <- which.max(logl)
  ww <- wl[[winner]]
  ss <- GG[, winner]
  names(ss) <- paste("shift", 1:(ng - 1), sep = "")
  ww$parameters <- append(ww$parameters, ss)
  ww$all.logl <- logl
  ww$GG <- GG
  out<-ww
    }

  return(out)

}
