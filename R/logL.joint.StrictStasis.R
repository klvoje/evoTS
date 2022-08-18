# Function from paleoTS (License:	GPL-3)

logL.joint.StrictStasis<- function (p, y)
{
  theta<- p[1]
  n<- length(y$mm)
  VV<- diag(y$vv/y$nn)
  detV<- det(VV)
  invV<- solve(VV)
  M<- rep(theta, n)
  #S<- dmvnorm(x$mm, mean = M, sigma = VV, log = TRUE)
  S<- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VV, log = TRUE)

  return(S)
}
