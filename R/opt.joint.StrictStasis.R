# Function from paleoTS (License:	GPL-3)

opt.joint.StrictStasis<- function (y, pool = TRUE, hess = FALSE)
{
  if (pool)
    y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)

  cl = list(fnscale = -1)

  p0 <- mean(y$mm)
  names(p0)<- "theta"
  w <- optim(p0, fn = logL.joint.StrictStasis, control = cl, method = "Brent", lower=min(y$mm), upper=max(y$mm),
             hessian = hess, y = y)
  names(w$par)<- "theta"
  if (hess)
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  wc<- paleoTS::as.paleoTSfit(logL=w$value, parameters=w$par, modelName="StrictStasis", method="Joint", K=1, n=length(y$mm), se=w$se)
  return(wc)
}
