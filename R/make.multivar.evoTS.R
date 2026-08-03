#'@title Makes a multivariate data set of
#'
#'@description Function to make a multivariate data set consisting of two or more evolutionary sequences (time-series).
#'
#'@param ... two or more univariate evolutionary sequences (time-series) in the format used by paleoTS, passed as individual arguments.
#'@details See the function as.paleoTS for details. See also read.paleoTS, which is often a more convenient way for getting the relevant data from text files.
#'@return a multivariate evoTS object that can be analysed with functions fitting multivariate models (e.g. fit.multivariate.OU, fit.multivariate.URW)
#'
#'@author Kjetil Lysne Voje
#'
#'@export
#
#'@examples
#'## Generate two evolutionary sequences (time-series)
#'x1 <- paleoTS::sim.GRW(60)
#'x2 <- paleoTS::sim.GRW(60)
#'
#'
#'## Make a multivariate data set
#'x1_x2<-make.multivar.evoTS(x1, x2)
#'

make.multivar.evoTS<-function (...){

  all.data <- list(...)

  xx <- do.call(cbind, lapply(all.data, `[[`, "mm"))
  vv <- do.call(cbind, lapply(all.data, `[[`, "vv"))
  nn <- do.call(cbind, lapply(all.data, `[[`, "nn"))
  tt <- do.call(cbind, lapply(all.data, `[[`, "tt"))

  y <- list(xx = xx, vv = vv, nn = nn, tt = tt)

  return(y)
}
