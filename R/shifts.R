# Function from paleoTS (License:	GPL-3)

shifts<-function (ns, ng, minb = 5) 
{
  aa <- combn(ns, ng - 1)
  if (ng == 2) 
    ok <- aa > minb & aa <= ns - minb + 1
  if (ng > 2) {
    daa <- apply(aa, 2, diff)
    if (ng > 3) 
      mdaa <- apply(daa, 2, min)
    else mdaa <- daa
    ok1 <- mdaa >= minb
    ok2 <- aa[1, ] > minb
    ok3 <- aa[ng - 1, ] <= ns - minb + 1
    ok <- ok1 & ok2 & ok3
  }
  ret <- aa[, ok]
  if (ng == 2) 
    ret <- matrix(ret, nrow = 1)
  if (is.null(dim(ret))) 
    ret <- matrix(ret, ncol = 1)
  if (sum(ok) > 0) 
    return(ret)
  else return(NULL)
}