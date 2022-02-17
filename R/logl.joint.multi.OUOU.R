#' @title Log-likelihoods for evolutionary models
#'
#' @description Returns log-likelihood for a multivariate Ornstein-Uhlenbeck model.
#'
#' @param init.par initial (starting) parameters values
#'
#' @param yy a multivariate evoTS object
#'
#' @param A.matrix the pull matrix.
#'
#' @param R.matrix the drift matrix..
#'
#' @details In general, users will not be access these functions directly, but instead use the optimization functions, which use these functions to find the best-supported parameter values.
#'
#'@return The log-likelihood of the parameter estimates, given the data.
#'
#'@author Kjetil Lysne Voje


logL.joint.multi.OUOU<-function (init.par, yy, A.matrix, R.matrix){

  m <-ncol(yy$xx) # number of traits
  X <- yy$xx # Character matrix with dimensions n * m
  y <- as.matrix(as.vector(X)) # Vectorized version of X

  if(A.matrix=="diag" & R.matrix=="diag"){

    ### Eigenvalue decomposition of A ###
    A<-diag(c(init.par[1:m]))
    P<-eigen(A)$vectors
    D<-diag(eigen(A)$values)

  ### The R (drift) matrix ###
  Chol<-diag(init.par[(m+1):(m*2)])

  ### Theta (optimal trait values) ###
  optima<-c(init.par[((m*2)+1):(m*3)])

  ### The ancestral trait values ###
  anc<-c(init.par[((m*3)+1):(m*4)])
  }


  if(A.matrix=="diag" & R.matrix=="symmetric"){

    ### Eigenvalue decomposition of A ###
    A<-diag(c(init.par[1:m]))
    P<-eigen(A)$vectors
    D<-diag(eigen(A)$values)

    ### The R (drift) matrix ###
    Chol<-diag(init.par[(m+1):(m*2)])
    nr.off.diag<-upper.tri(Chol)
    l.upp.tri<-length(nr.off.diag[nr.off.diag==TRUE])
    Chol[upper.tri(Chol)] <- init.par[((length(diag(A))*2)+1):((length(diag(A))*2)+l.upp.tri)]

    ### Theta (optimal trait values) ###
    optima<-c(init.par[((length(diag(A))*2)+l.upp.tri+1):((length(diag(A))*2)+l.upp.tri+m)])

    ### The ancestral trait values ###
    anc<-c(init.par[((length(diag(A))*2)+l.upp.tri+1+m):((length(diag(A))*2)+l.upp.tri+m+m)])
  }

  if(A.matrix=="full" & R.matrix=="diag"){

    ### Eigenvalue decomposition of A ###
    A<-diag(c(init.par[1:m]))
    nr.off.diag<-upper.tri(A)
    l.upp.tri<-length(nr.off.diag[nr.off.diag==TRUE])
    A[upper.tri(A)] <- init.par[(length(diag(A))+1):(length(diag(A))+l.upp.tri)]
    A[lower.tri(A)] <- init.par[(length(diag(A))+l.upp.tri+1):(length(diag(A))+l.upp.tri+l.upp.tri)]
    P<-eigen(A)$vectors
    D<-diag(eigen(A)$values)

    ### The R (drift) matrix ###
    Chol<-diag(init.par[(length(diag(A))+l.upp.tri+l.upp.tri+1):(length(diag(A))+l.upp.tri+l.upp.tri+(length(diag(A))))])

    ### Theta (optimal trait values) ###
    optima<-c(init.par[(length(diag(A))+l.upp.tri+l.upp.tri+(length(diag(A)))+1):(length(diag(A))+l.upp.tri+l.upp.tri+(length(diag(A)))+length(diag(A)))])

    ### The ancestral trait values ###
    anc<-c(init.par[(length(diag(A))+l.upp.tri+l.upp.tri+(length(diag(A)))+length(diag(A))+1):(length(diag(A))+l.upp.tri+l.upp.tri+(length(diag(A)))+length(diag(A))++length(diag(A)))])

  }

  if(A.matrix=="full" & R.matrix=="symmetric"){

    ### Eigenvalue decomposition of A ###
    A<-diag(c(init.par[1:m]))
    nr.off.diag<-upper.tri(A)
    l.upp.tri<-length(nr.off.diag[nr.off.diag==TRUE])
    A[upper.tri(A)] <- init.par[(length(diag(A))+1):(length(diag(A))+l.upp.tri)]
    A[lower.tri(A)] <- init.par[(length(diag(A))+l.upp.tri+1):(length(diag(A))+l.upp.tri+l.upp.tri)]
    P<-eigen(A)$vectors
    D<-diag(eigen(A)$values)

    ### The R (drift) matrix ###
    Chol<-diag(init.par[(length(diag(A))+l.upp.tri+l.upp.tri+1):(length(diag(A))+l.upp.tri+l.upp.tri+length(diag(A)))])
    nr.off.diag.R<-upper.tri(Chol)
    l.upp.tri.R<-length(nr.off.diag.R[nr.off.diag.R==TRUE])
    Chol[upper.tri(Chol)] <- init.par[(length(diag(A))+l.upp.tri+l.upp.tri+length(diag(A))+1):(length(diag(A))+l.upp.tri+l.upp.tri+length(diag(A))+l.upp.tri.R)]

    ### Theta (optimal trait values) ###
    optima<-c(init.par[(length(diag(A))+l.upp.tri+l.upp.tri+length(diag(A))+l.upp.tri.R+1):(length(diag(A))+l.upp.tri+l.upp.tri+length(diag(A))+l.upp.tri.R+length(diag(A)))])

    ### The ancestral trait values ###
    anc<-c(init.par[(length(diag(A))+l.upp.tri+l.upp.tri+length(diag(A))+l.upp.tri.R+length(diag(A))+1):(length(diag(A))+l.upp.tri+l.upp.tri+length(diag(A))+l.upp.tri.R+length(diag(A))+length(diag(A)))])

  }

  if(A.matrix=="upper.tri" & R.matrix=="diag"){

    ### Eigenvalue decomposition of A ###
    A<-diag(c(init.par[1:m]))
    nr.off.diag<-upper.tri(A)
    l.upp.tri<-length(nr.off.diag[nr.off.diag==TRUE])
    A[upper.tri(A)] <- init.par[(length(diag(A))+1):(length(diag(A))+l.upp.tri)]
    P<-eigen(A)$vectors
    D<-diag(eigen(A)$values)

    ### The R (drift) matrix ###
    Chol<-diag(init.par[(length(diag(A))+l.upp.tri+1):((length(diag(A))+l.upp.tri+m))])

    ### Theta (optimal trait values) ###
    optima<-c(init.par[(length(diag(A))+l.upp.tri+m+1):((length(diag(A))+l.upp.tri+m+m))])

    ### The ancestral trait values ###
    anc<-c(init.par[(length(diag(A))+l.upp.tri+m+m+1):(length(diag(A))+l.upp.tri+m+m+m)])

  }

  if(A.matrix=="upper.tri" & R.matrix=="symmetric"){

    ### Eigenvalue decomposition of A ###
    A<-diag(c(init.par[1:m]))
    nr.off.diag.A<-upper.tri(A)
    l.upp.tri.A<-length(nr.off.diag.A[nr.off.diag.A==TRUE])
    A[upper.tri(A)] <- init.par[(length(diag(A))+1):(length(diag(A))+l.upp.tri.A)]
    P<-eigen(A)$vectors
    D<-diag(eigen(A)$values)

    ### The R (drift) matrix ###
    Chol<-diag(init.par[(length(diag(A))+l.upp.tri.A+1):(length(diag(A))+l.upp.tri.A+m)])
    nr.off.diag.R<-upper.tri(Chol)
    l.upp.tri.R<-length(nr.off.diag.R[nr.off.diag.R==TRUE])
    Chol[upper.tri(Chol)] <- init.par[(length(diag(A))+l.upp.tri.A+m+1):(length(diag(A))+l.upp.tri.A+m+l.upp.tri.R)]

    ### Theta (optimal trait values) ###
    optima<-c(init.par[(length(diag(A))+l.upp.tri.A+m+l.upp.tri.R+1):(length(diag(A))+l.upp.tri.A+m+l.upp.tri.R+m)])

    ### The ancestral trait values ###
    anc<-c(init.par[(length(diag(A))+l.upp.tri.A+m+l.upp.tri.R+m+1):(length(diag(A))+l.upp.tri.A+m+l.upp.tri.R+m+m)])
  }

  if(A.matrix=="lower.tri" & R.matrix=="diag"){

    ### Eigenvalue decomposition of A ###
    A<-diag(c(init.par[1:m]))
    nr.off.diag<-lower.tri(A)
    l.low.tri<-length(nr.off.diag[nr.off.diag==TRUE])
    A[lower.tri(A)] <- init.par[(length(diag(A))+1):(length(diag(A))+l.low.tri)]
    P<-eigen(A)$vectors
    D<-diag(eigen(A)$values)

    ### The R (drift) matrix ###
    Chol<-diag(init.par[(length(diag(A))+l.low.tri+1):((length(diag(A))+l.low.tri+m))])

    ### Theta (optimal trait values) ###
    optima<-c(init.par[(length(diag(A))+l.low.tri+m+1):((length(diag(A))+l.low.tri+m+m))])

    ### The ancestral trait values ###
    anc<-c(init.par[(length(diag(A))+l.low.tri+m+m+1):(length(diag(A))+l.low.tri+m+m+m)])

  }

  if(A.matrix=="lower.tri" & R.matrix=="symmetric"){

    ### Eigenvalue decomposition of A ###
    A<-diag(c(init.par[1:m]))
    nr.off.diag.A<-lower.tri(A)
    l.low.tri.A<-length(nr.off.diag.A[nr.off.diag.A==TRUE])
    A[lower.tri(A)] <- init.par[(length(diag(A))+1):(length(diag(A))+l.low.tri.A)]
    P<-eigen(A)$vectors
    D<-diag(eigen(A)$values)

    ### The R (drift) matrix ###
    Chol<-diag(init.par[(length(diag(A))+l.low.tri.A+1):(length(diag(A))+l.low.tri.A+m)])
    nr.off.diag.R<-upper.tri(Chol)
    l.upp.tri.R<-length(nr.off.diag.R[nr.off.diag.R==TRUE])
    Chol[upper.tri(Chol)] <- init.par[(length(diag(A))+l.low.tri.A+m+1):(length(diag(A))+l.low.tri.A+m+l.upp.tri.R)]

    ### Theta (optimal trait values) ###
    optima<-c(init.par[(length(diag(A))+l.low.tri.A+m+l.upp.tri.R+1):(length(diag(A))+l.low.tri.A+m+l.upp.tri.R+m)])

    ### The ancestral trait values ###
    anc<-c(init.par[(length(diag(A))+l.low.tri.A+m+l.upp.tri.R+m+1):(length(diag(A))+l.low.tri.A+m+l.upp.tri.R+m+m)])
  }

  if(A.matrix=="OUBM"){

    ### Eigenvalue decomposition of A ###
    A<-diag(c(init.par[1:(m-1)],0.000001))
    nr.off.diag.A<-upper.tri(A)
    l.upp.tri.A<-length(nr.off.diag.A[nr.off.diag.A==TRUE])
    A[upper.tri(A)] <- init.par[(m):(m-1+l.upp.tri.A)]
    P<-eigen(A)$vectors
    D<-diag(eigen(A)$values)

    ### The R (drift) matrix ###
    Chol<-diag(init.par[(m+l.upp.tri.A):(m+m+l.upp.tri.A-1)])

    ### Theta (optimal trait values) ###
    optima<-c(init.par[(m+m+l.upp.tri.A):(m+m+m+l.upp.tri.A-1)])

    ### The ancestral trait values ###
    anc<-c(init.par[(m+m+m+l.upp.tri.A):(m+m+m+m+l.upp.tri.A-1)])
  }

  ### Define and rescale the time vector to unit length ###
   time<-yy$tt[,1]/max(yy$tt[,1])

  # Make a time matrix
   tmp.matrix<-matrix(0,ncol=length(yy$vv[,1]), nrow=length(yy$vv[1,]))
  for (i in 1:ncol(A)){
    tmp.matrix[i,]<-time
  }
  time<-tmp.matrix

  ### Calculate exponent of eigenvalues from the A decomposition ###

  exp_eigenvalues<-array(data=NA, dim=c(m, m, length(time[1,])))

  for (i in 1:length(time[1,])){
    exp_eigenvalues_tmp<-rep(NA, m)
      for (j in 1: m){
        exp_eigenvalues_tmp[j]<-exp(-diag(D)[j]*time[1,i])
      }
      exp_eigenvalues[,,i]<-diag(exp_eigenvalues_tmp)
  }


  ### Calculate the expected trait evolution given A, anc, theta and time. ###
  if(A.matrix=="full" || A.matrix=="upper.tri" || A.matrix=="lower.tri" || A.matrix=="OUBM")
    {
    M<-matrix(NA, ncol = length(time[1,]), nrow= m)
    for (i in 1:length(time[1,]))
      {
      M[,i]<-((P%*%exp_eigenvalues[,,i]%*%solve(P))%*%anc) + (diag(c(rep(1,m)))- (P%*%exp_eigenvalues[,,i]%*%solve(P)))%*%optima
      }
    }

  if(A.matrix=="diag")
    {
  M<-exp(-(P%*%D%*%solve(P))%*%time)*anc + optima*(1- exp(-(P%*%D%*%solve(P))%*%time))
    }
   M<-c(t(M)) # vectorize M


  ### Compute the integral of the variance-covariance matrix (eq. 8 and 9 (page 3) from Suppl. from Clavel et al. 2015) ###

  tmp.VV<-array(data=NA, dim=c(ncol=length(time[1,]), nrow=length(time[1,]), (ncol(A)*ncol(A)))) # Make a list that contains the block matrices in the VCOV.

  # Create time matrices (ta and tij)
  ff <- function(a, b) abs(a - b)
  tij<-outer(as.vector(time[1,]), as.vector(time[1,]), ff) # tij -> time from species j to the most common ancestor of species i and j.
  ta<-outer(as.vector(time[1,]), as.vector(time[1,]),pmin) #Ta -> time from the first sample to the most recent common ancestor of i and j.

  right.side<-solve(P)%*% (Chol %*% t(Chol)) %*% (t(solve(P))) # the right side of the expression relative to the Hadamard product

  left.side<-matrix(NA, nrow=nrow(A), ncol=ncol(A)) # the left side of the expression relative to the Hadamard product
  for (i in 1:length(time[1,])){
    for (j in 1:length(time[1,])){
      for (k in 1:ncol(A)){
        for (l in 1:ncol(A)){
          left.side[k,l]<-(1-exp(-(diag(D)[k]+diag(D)[l])*ta[i,j]))*(1/(diag(D)[k]+diag(D)[l]))
          }
      }
     # print(left.side)
      left.right<-left.side*right.side # Hadamard product between left and right side
      integ<-P%*%left.right%*%t(P) # The whole integral (except the matrix exponential)

        exp_eigenvalues_2<-rep(NA,  m)
        for (m in 1:m){
        exp_eigenvalues_2[m]<-exp(-diag(D)[m]*tij[i,j])
        }

          tmp<-integ%*%(t(P%*%(diag(c(exp_eigenvalues_2)))%*%solve(P))) # The whole integral (including the matrix exponential)
          vector.tmp<-as.vector(tmp) # vectorization of the trait*trait matrix
          for (k in 1:(ncol(A)*ncol(A))){
          tmp.VV[i,j,k]<-vector.tmp[k] # place the elements of the integral in each block matrix.
        }
      }
  }

  ### Compile the varcovar (VV3) matrix from the list of block matrices (tmp.VV) ###

  VV3<-matrix(0, ncol=(ncol(tmp.VV[,,1])*ncol(A)), nrow=(ncol(tmp.VV[,,1])*ncol(A))) # Make an empty VCOV matrix.

  List<-list()
  for(i in 1:length(tmp.VV[1,1,]))
  {
      List[[i]] <- tmp.VV[,,i] # Make each block matrix a separate element in the list "List"
  }

  #List[[3]]<-t( List[[3]])

from.boundary<-seq(1,ncol(VV3), ncol(VV3)/ncol(A)) # create vectors defining the start...
to.boundary<-(from.boundary-1)+length(tmp.VV[,1,1]) # and end of where in the VCOV matrix block matrices in List should be placed.
from<-seq(1,length(tmp.VV[1,1,]), length(tmp.VV[1,1,])/ncol(A)) # create vectors defining the start and and end of which list in List that should be placed in VCOV.
to<-(from-1)+ncol(A)

for (i in 1:ncol(A)){
  VV3[(from.boundary[i]:to.boundary[i]),]<-do.call(cbind, List[from[i]: to[i]]) # Make the VCOV matrix by binding together block matrices from List
}

 #print(VV3)
#### Add estimation error to the diagonal ###

  tmp.matrix<-matrix(0,ncol=length(yy$vv[,1]), nrow=length(yy$vv[1,])) # Make empty matrix for sampling error
  for (i in 1:ncol(A)){
    tmp.matrix[i,]<-c(yy$vv[,i]/yy$nn[,i]) # add sampling error (sample variance divided by sample size)
  }

  diag(VV3) <- diag(VV3) + as.vector(t(tmp.matrix)) # Add sampling error to the diagonal of VCOV
  VV3<-(VV3+t(VV3))/2 #Make sure round-off errors are not present in VV3

  S <- mvtnorm::dmvnorm(t(y), mean = M, sigma = VV3, log = TRUE)

}
