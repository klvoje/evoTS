shift2gg<-function (ss, ns) 
{
  z <- c(0, ss, ns + 1)
  cc <- cut(1:ns, breaks = z, right = FALSE)
  gg <- as.numeric(cc)
  return(gg)
}