#' @title Log-likelihoods for evolutionary models
#'
#' @description Returns log-likelihood for a multivariate Ornstein-Uhlenbeck model with used defined A and R matrices..
#'
#' @param init.par initial (starting) parameters values
#'
#' @param yy a multivariate evoTS object
#'
#' @param A.user the pull matrix.
#'
#' @param R.user the drift matrix.
#'
#' @param locations.A location (row and column) of parameters (elements) in the A matrix that is estimated
#'
#' @param location.diag.A location (row and column) of parameters (elements) in the diagonal of the A matrix that is estimated
#'
#' @param location.upper.tri.A location (row and column) of parameters (elements) in the upper triangle of the A matrix that is estimated
#'
#' @param location.lower.tri.A location (row and column) of parameters (elements) in the lower triangle of the A matrix that is estimated
#'
#' @param locations.R location (row and column) of parameters (elements) in the R matrix that is estimated
#'
#' @param location.diag.R location (row and column) of parameters (elements) in the diagonal of the R matrix that is estimated
#'
#' @param location.upper.tri.R location (row and column) of parameters (elements) in the upper triangle of the R matrix that is estimated
#'
#' @details In general, users will not be access these functions directly, but instead use the optimization functions, which use these functions to find the best-supported parameter values.
#'
#'@return The log-likelihood of the parameter estimates, given the data.
#'
#'@author Kjetil Lysne Voje


logL.joint.multi.OUOU.user<-function (init.par, yy, A.user, R.user, locations.A, location.diag.A, location.upper.tri.A, location.lower.tri.A,
                                      locations.R, location.diag.R, location.upper.tri.R){

  m <-ncol(yy$xx) # number of traits
  X <- yy$xx # Character matrix with dimensions n * m
  y <- as.matrix(as.vector(X)) # Vectorized version of X

  A<-diag(rep(0.00000000001,m))
  #A[locations.A[location.diag.A],locations.A[location.diag.A]]<- diag(c(init.par[1:length(location.diag.A)]))
  for (i in 1:length(location.diag.A)){
  A[locations.A[location.diag.A][i],locations.A[location.diag.A][i]]<- init.par[i]
  }
  
  if (pracma::isempty(location.upper.tri.A)==FALSE)
    {
    
    for (i in 1:length(location.upper.tri.A)){
      A[locations.A[,1][location.upper.tri.A][i],locations.A[,2][location.upper.tri.A][i]]<- init.par[(length(location.diag.A)+i)]
    }
  } else location.upper.tri.A<-NULL
  
  #  A[locations.A[,1][location.upper.tri.A],locations.A[,2][location.upper.tri.A]]<-init.par[(length(location.diag.A)+1):(length(location.diag.A)+length(location.upper.tri.A))]
   # } else location.upper.tri.A<-NULL

  if (pracma::isempty(location.lower.tri.A)==FALSE)
    {
    for (i in 1:length(location.lower.tri.A)){
  A[locations.A[,1][location.lower.tri.A][i],locations.A[,2][location.lower.tri.A][i]]<-init.par[(length(location.diag.A)+length(location.upper.tri.A)+i)]
  } 
    }else location.lower.tri.A<-NULL

  P<-eigen(A)$vectors
  D<-diag(eigen(A)$values)


  Chol<-diag(rep(0,m))
  for (i in 1:length(location.diag.R)){
  Chol[locations.R[location.diag.R][i],locations.R[location.diag.R][i]]<- init.par[(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+i)]
  }
  
  if (pracma::isempty(location.upper.tri.R)==FALSE)
  {
    for (i in 1:length(location.upper.tri.R)){
  Chol[locations.R[,1][location.upper.tri.R][i],locations.R[,2][location.upper.tri.R][i]]<-init.par[(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+length(location.diag.R)+i)]
    }
  } else location.upper.tri.R<-NULL

  ### Theta (optimal trait values) ###
  optima<-c(init.par[(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+length(location.diag.R)+length(location.upper.tri.R)+1):(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+length(location.diag.R)+length(location.upper.tri.R)+m)])

  ### The ancestral trait values ###
  anc<-c(init.par[(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+length(location.diag.R)+length(location.upper.tri.R)+m+1):(length(location.diag.A)+length(location.upper.tri.A)+length(location.lower.tri.A)+length(location.diag.R)+length(location.upper.tri.R)+m+m)])


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

    M<-matrix(NA, ncol = length(time[1,]), nrow= m)
    for (i in 1:length(time[1,]))
      {
      M[,i]<-((P%*%exp_eigenvalues[,,i]%*%solve(P))%*%anc) + (diag(c(rep(1,m)))- (P%*%exp_eigenvalues[,,i]%*%solve(P)))%*%optima
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

#### Add estimation error to the diagonal ###

  tmp.matrix<-matrix(0,ncol=length(yy$vv[,1]), nrow=length(yy$vv[1,])) # Make empty matrix for sampling error
  for (i in 1:ncol(A)){
    tmp.matrix[i,]<-c(yy$vv[,i]/yy$nn[,i]) # add sampling error (sample variance divided by sample size)
  }

  diag(VV3) <- diag(VV3) + as.vector(t(tmp.matrix)) # Add sampling error to the diagonal of VCOV
  VV3<-(VV3+t(VV3))/2 #Make sure round-off errors are not present in VV3

  S <- mvtnorm::dmvnorm(t(y), mean = M, sigma = VV3, log = TRUE)

}
