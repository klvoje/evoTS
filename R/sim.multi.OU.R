#' @title Simulate multivariate Ornstein-Uhlenbeck evolutionary sequence data sets
#'
#' @description Function to simulate a multivariate Ornstein-Uhlenbeck evolutionary sequence data set.
#'
#' @param ns number of samples in time-series
#'
#' @param anc the ancestral trait values
#'
#' @param optima the optimal trait values
#'
#' @param A the pull matrix.
#'
#' @param R the drift matrix
#'
#' @param vp within-population trait variance
#'
#' @param nn 	vector of the number of individuals in each sample (identical sample sizes for all time-series is assumed)
#'
#' @param tt 	vector of sample ages, increases from oldest to youngest
#'
#'@return A multivariate evolutionary sequence (time-series) data set.
#'
#'@note The Ornstein Uhlenbeck model is reduced to an Unbiased Random Walk when the alpha parameter is zero. It is therefore possible to let a trait evolve as an Unbiased Random Walk by setting the diagonal element for that trait to a value close to zero (e.g. 1e-07). Elements in the diagonal of A cannot be exactly zero as this will result in a singular variance-covariance matrix.
#'
#'@author Kjetil Lysne Voje
#'
#'@export
#'
#'@examples
#'##Define the A and R matrices
#'
#'A_matrix<-matrix(c(4,-2,0,3), nrow=2, byrow = TRUE)
#'R_matrix<-matrix(c(4,0.2,0.2,4), nrow=2, byrow = TRUE)
#'
#'## Generate an evoTS object by simulating a multivariate dataset
#'data_set<-sim.multi.OU(40, optima = c(1.5,2),A=A_matrix , R = R_matrix)
#'
#'## plot the data
#'plotevoTS.multivariate(data_set)
#'
sim.multi.OU<-function(ns = 30, anc = c(0,0), optima = c(3, 2),
                       A= matrix(c(7,0,0,6), nrow=2, byrow = TRUE), R = matrix(c(0.05,0,0,0.05), nrow=2, byrow = TRUE),
                       vp = 0.1, nn = rep(30, ns), tt = 0:(ns - 1)){
  m<-ncol(A)

  A.matrix<-A
  #A.matrix[1,2]<-A[1,2]
  P<-eigen(A.matrix)$vectors
  D<-diag(eigen(A.matrix)$values)
  Chol<-chol(R)

  # Make a time matrix
  time<-matrix(nrow = ncol(A), ncol = ns)
  for (i in 1:ncol(A)){
    time[i,]<-tt/max(tt)
  }

  tmp.VV<-array(data=NA, dim=c(ncol=length(time[1,]), nrow=length(time[1,]), (ncol(A)*ncol(A)))) # Make a list that contains the block matrices in the VCOV.

  ff <- function(a, b) abs(a - b)
  tij<-outer(as.vector(time[1,]), as.vector(time[1,]), ff) # tij -> time from species j to the most common ancestor of species i and j.
  #tij<-tij[1,]
  ta<-outer(as.vector(time[1,]), as.vector(time[1,]),pmin) #Ta -> time from the first sample to the most recent common ancestor of i and j.
#  ta<-diag(ta)
  MM <- matrix(nrow = ncol(A), ncol = ns)
  mm <- matrix(nrow = ncol(A), ncol = ns)
  vv <- matrix(nrow = ncol(A), ncol = ns)
  MM[,1] <- anc

  for (i in 1:m){
    x <- rnorm(nn[i], mean = MM[i], sd = sqrt(vp))
    mm[i] <- mean(x)
    vv[i] <- var(x)
  }

  traits<-matrix(NA, ncol = length(time[1,]), nrow= m)
  for (i in 1:length(time[1,])){
    traits[,i]<-((P%*%diag(c(exp(-diag(D)[1]*time[1,i]),exp(-diag(D)[2]*time[1,i])))%*%solve(P))%*%anc) + (matrix(c(1,0,0,1), ncol=2, byrow=2)- (P%*%diag(c(exp(-diag(D)[1]*time[1,i]),exp(-diag(D)[2]*time[1,i])))%*%solve(P)))%*%optima

  }

  varcovar<-array(data =NA, dim=c(m,m,ns))

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

  VV3<-matrix(0, ncol=(ncol(tmp.VV[,,1])*ncol(A)), nrow=(ncol(tmp.VV[,,1])*ncol(A))) # Make an empty VCOV matrix.

  List<-list()
  for(i in 1:length(tmp.VV[1,1,]))
  {
    List[[i]] <- tmp.VV[,,i] # Make each block matrix a separate element in the list "List"
  }

  from.boundary<-seq(1,ncol(VV3), ncol(VV3)/ncol(A)) # create vectors defining the start...
  to.boundary<-(from.boundary-1)+length(tmp.VV[,1,1]) # and end of where in the VCOV matrix block matrices in List should be placed.
  from<-seq(1,length(tmp.VV[1,1,]), length(tmp.VV[1,1,])/ncol(A)) # create vectors defining the start and and end of which list in List that should be placed in VCOV.
  to<-(from-1)+ncol(A)

  for (i in 1:ncol(A)){
    VV3[(from.boundary[i]:to.boundary[i]),]<-do.call(cbind, List[from[i]: to[i]]) # Make the VCOV matrix by binding together block matrices from List
  }

  traits_x<-c(t(traits))

      MM<-MASS::mvrnorm(n = 1, mu=traits_x, Sigma=VV3)#, tol = 1e-6, empirical = TRUE, EISPACK = FALSE)
      MM<-matrix(MM, m,ns, byrow=TRUE)
      for (i in 2:ns) {
      for (j in 1:m){

         x <- MASS::mvrnorm(nn[j], mu = MM[j,i], Sigma = sqrt(vp)) #NEW. rnorm can't handle irrational numbers.
         #x <- ?rnorm(nn[j], mean = MM[j,i], sd = sqrt(vp))
      mm[j,i] <- mean(as.numeric(x)) #as.numeric is needed in case irrational numbers are present in x.
      vv[j,i] <- var(x)
    }

        
  }


  List<-list()
  for (i in 1:m){
    List[[i]]<-paleoTS::as.paleoTS(mm = mm[i,], vv = vv[i,], nn = nn, tt = time[i,], MM = MM[i,], label = "Created by sim.multi.OU()", reset.time = FALSE)
  }

  if (m==2) yy<-make.multivar.evoTS(List[[1]], List[[2]])
  if (m==3) yy<-make.multivar.evoTS(List[[1]], List[[2]], List[[3]])
  if (m==4) yy<-make.multivar.evoTS(List[[1]], List[[2]], List[[3]], List[[4]])
  if (m==5) yy<-make.multivar.evoTS(List[[1]], List[[2]], List[[3]], List[[4]], List[[5]])


  return(yy)
}
