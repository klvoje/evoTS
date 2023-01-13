#'@title Makes a multivariate data set of
#'
#'@description Function to make a multivariate data set consisting of two or more evolutionary sequences (time-series).
#'
#'@param evoTS.1 an univariate evolutionary sequences (time-series) on the format used in paleoTS
#'@param evoTS.2 an univariate evolutionary sequences (time-series) on the format used in paleoTS
#'@param evoTS.3 an univariate evolutionary sequences (time-series) on the format used in paleoTS (optional)
#'@param evoTS.4 an univariate evolutionary sequences (time-series) on the format used in paleoTS (optional)
#'@param evoTS.5 an univariate evolutionary sequences (time-series) on the format used in paleoTS (optional)
#'@param evoTS.6 an univariate evolutionary sequences (time-series) on the format used in paleoTS (optional)
#'@param evoTS.7 an univariate evolutionary sequences (time-series) on the format used in paleoTS (optional)
#'@param evoTS.8 an univariate evolutionary sequences (time-series) on the format used in paleoTS (optional)
#'@param evoTS.9 an univariate evolutionary sequences (time-series) on the format used in paleoTS (optional)
#'@param evoTS.10 an univariate evolutionary sequences (time-series) on the format used in paleoTS (optional)
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

make.multivar.evoTS<-function (evoTS.1=NULL, evoTS.2=NULL, evoTS.3=NULL, evoTS.4=NULL, evoTS.5=NULL, evoTS.6=NULL, evoTS.7=NULL, evoTS.8=NULL, evoTS.9=NULL, evoTS.10=NULL){

  if (is.numeric(evoTS.3$mm[1]) == TRUE) all.data<-list(evoTS.1, evoTS.2, evoTS.3) else all.data<-list(evoTS.1, evoTS.2)
  if (is.numeric(evoTS.4$mm[1]) == TRUE) all.data<-list(evoTS.1, evoTS.2, evoTS.3, evoTS.4) else all.data<-list(evoTS.1, evoTS.2, evoTS.3)
  if (is.numeric(evoTS.5$mm[1]) == TRUE) all.data<-list(evoTS.1, evoTS.2, evoTS.3, evoTS.4, evoTS.5) else all.data<-list(evoTS.1, evoTS.2, evoTS.3, evoTS.4)
  if (is.numeric(evoTS.6$mm[1]) == TRUE) all.data<-list(evoTS.1, evoTS.2, evoTS.3, evoTS.4, evoTS.5, evoTS.6) else all.data<-list(evoTS.1, evoTS.2, evoTS.3, evoTS.4, evoTS.5)
  if (is.numeric(evoTS.7$mm[1]) == TRUE) all.data<-list(evoTS.1, evoTS.2, evoTS.3, evoTS.4, evoTS.5, evoTS.6, evoTS.7) else all.data<-list(evoTS.1, evoTS.2, evoTS.3, evoTS.4, evoTS.5, evoTS.6)
  if (is.numeric(evoTS.8$mm[1]) == TRUE) all.data<-list(evoTS.1, evoTS.2, evoTS.3, evoTS.4, evoTS.5, evoTS.6, evoTS.7, evoTS.8) else all.data<-list(evoTS.1, evoTS.2, evoTS.3, evoTS.4, evoTS.5, evoTS.6, evoTS.7)
  if (is.numeric(evoTS.9$mm[1]) == TRUE) all.data<-list(evoTS.1, evoTS.2, evoTS.3, evoTS.4, evoTS.5, evoTS.6, evoTS.7, evoTS.8, evoTS.9) else all.data<-list(evoTS.1, evoTS.2, evoTS.3, evoTS.4, evoTS.5, evoTS.6, evoTS.7, evoTS.8)
  if (is.numeric(evoTS.10$mm[1]) == TRUE) all.data<-list(evoTS.1, evoTS.2, evoTS.3, evoTS.4, evoTS.5, evoTS.6, evoTS.7, evoTS.8, evoTS.9, evoTS.10) else all.data<-list(evoTS.1, evoTS.2, evoTS.3, evoTS.4, evoTS.5, evoTS.6, evoTS.7, evoTS.8, evoTS.9)

  number.of.traits<-length(all.data)

  if(number.of.traits ==2)  {
    xx<-cbind(all.data[[1]]$mm, all.data[[2]]$mm)
    vv<-cbind(all.data[[1]]$vv, all.data[[2]]$vv)
    nn<-cbind(all.data[[1]]$nn, all.data[[2]]$nn)
    tt<-cbind(all.data[[1]]$tt, all.data[[2]]$tt)
  }

  if(number.of.traits ==3)  {
    xx<-cbind(all.data[[1]]$mm, all.data[[2]]$mm, all.data[[3]]$mm)
    vv<-cbind(all.data[[1]]$vv, all.data[[2]]$vv, all.data[[3]]$vv)
    nn<-cbind(all.data[[1]]$nn, all.data[[2]]$nn, all.data[[3]]$nn)
    tt<-cbind(all.data[[1]]$tt, all.data[[2]]$tt, all.data[[3]]$tt)
  }

  if(number.of.traits ==4)  {
    xx<-cbind(all.data[[1]]$mm, all.data[[2]]$mm, all.data[[3]]$mm, all.data[[4]]$mm)
    vv<-cbind(all.data[[1]]$vv, all.data[[2]]$vv, all.data[[3]]$vv, all.data[[4]]$vv)
    nn<-cbind(all.data[[1]]$nn, all.data[[2]]$nn, all.data[[3]]$nn, all.data[[4]]$nn)
    tt<-cbind(all.data[[1]]$tt, all.data[[2]]$tt, all.data[[3]]$tt, all.data[[4]]$tt)
  }

  if(number.of.traits ==5)  {
    xx<-cbind(all.data[[1]]$mm, all.data[[2]]$mm, all.data[[3]]$mm, all.data[[4]]$mm, all.data[[5]]$mm)
    vv<-cbind(all.data[[1]]$vv, all.data[[2]]$vv, all.data[[3]]$vv, all.data[[4]]$vv, all.data[[5]]$vv)
    nn<-cbind(all.data[[1]]$nn, all.data[[2]]$nn, all.data[[3]]$nn, all.data[[4]]$nn, all.data[[5]]$nn)
    tt<-cbind(all.data[[1]]$tt, all.data[[2]]$tt, all.data[[3]]$tt, all.data[[4]]$tt, all.data[[5]]$tt)
  }

  if(number.of.traits ==6)  {
    xx<-cbind(all.data[[1]]$mm, all.data[[2]]$mm, all.data[[3]]$mm, all.data[[4]]$mm, all.data[[5]]$mm, all.data[[6]]$mm)
    vv<-cbind(all.data[[1]]$vv, all.data[[2]]$vv, all.data[[3]]$vv, all.data[[4]]$vv, all.data[[5]]$vv, all.data[[6]]$vv)
    nn<-cbind(all.data[[1]]$nn, all.data[[2]]$nn, all.data[[3]]$nn, all.data[[4]]$nn, all.data[[5]]$nn, all.data[[6]]$nn)
    tt<-cbind(all.data[[1]]$tt, all.data[[2]]$tt, all.data[[3]]$tt, all.data[[4]]$tt, all.data[[5]]$tt, all.data[[6]]$tt)
  }

  if(number.of.traits ==7)  {
    xx<-cbind(all.data[[1]]$mm, all.data[[2]]$mm, all.data[[3]]$mm, all.data[[4]]$mm, all.data[[5]]$mm, all.data[[6]]$mm, all.data[[7]]$mm)
    vv<-cbind(all.data[[1]]$vv, all.data[[2]]$vv, all.data[[3]]$vv, all.data[[4]]$vv, all.data[[5]]$vv, all.data[[6]]$vv, all.data[[7]]$vv)
    nn<-cbind(all.data[[1]]$nn, all.data[[2]]$nn, all.data[[3]]$nn, all.data[[4]]$nn, all.data[[5]]$nn, all.data[[6]]$nn, all.data[[7]]$nn)
    tt<-cbind(all.data[[1]]$tt, all.data[[2]]$tt, all.data[[3]]$tt, all.data[[4]]$tt, all.data[[5]]$tt, all.data[[6]]$tt, all.data[[7]]$tt)
  }

  if(number.of.traits ==8)  {
    xx<-cbind(all.data[[1]]$mm, all.data[[2]]$mm, all.data[[3]]$mm, all.data[[4]]$mm, all.data[[5]]$mm, all.data[[6]]$mm, all.data[[7]]$mm, all.data[[8]]$mm)
    vv<-cbind(all.data[[1]]$vv, all.data[[2]]$vv, all.data[[3]]$vv, all.data[[4]]$vv, all.data[[5]]$vv, all.data[[6]]$vv, all.data[[7]]$vv, all.data[[8]]$vv)
    nn<-cbind(all.data[[1]]$nn, all.data[[2]]$nn, all.data[[3]]$nn, all.data[[4]]$nn, all.data[[5]]$nn, all.data[[6]]$nn, all.data[[7]]$nn, all.data[[8]]$nn)
    tt<-cbind(all.data[[1]]$tt, all.data[[2]]$tt, all.data[[3]]$tt, all.data[[4]]$tt, all.data[[5]]$tt, all.data[[6]]$tt, all.data[[7]]$tt, all.data[[8]]$tt)
  }

  if(number.of.traits ==9)  {
    xx<-cbind(all.data[[1]]$mm, all.data[[2]]$mm, all.data[[3]]$mm, all.data[[4]]$mm, all.data[[5]]$mm, all.data[[6]]$mm, all.data[[7]]$mm, all.data[[8]]$mm, all.data[[9]]$mm)
    vv<-cbind(all.data[[1]]$vv, all.data[[2]]$vv, all.data[[3]]$vv, all.data[[4]]$vv, all.data[[5]]$vv, all.data[[6]]$vv, all.data[[7]]$vv, all.data[[8]]$vv, all.data[[9]]$vv)
    nn<-cbind(all.data[[1]]$nn, all.data[[2]]$nn, all.data[[3]]$nn, all.data[[4]]$nn, all.data[[5]]$nn, all.data[[6]]$nn, all.data[[7]]$nn, all.data[[8]]$nn, all.data[[9]]$nn)
    tt<-cbind(all.data[[1]]$tt, all.data[[2]]$tt, all.data[[3]]$tt, all.data[[4]]$tt, all.data[[5]]$tt, all.data[[6]]$tt, all.data[[7]]$tt, all.data[[8]]$tt, all.data[[9]]$tt)
  }

  if(number.of.traits ==10)  {
    xx<-cbind(all.data[[1]]$mm, all.data[[2]]$mm, all.data[[3]]$mm, all.data[[4]]$mm, all.data[[5]]$mm, all.data[[6]]$mm, all.data[[7]]$mm, all.data[[8]]$mm, all.data[[9]]$mm, all.data[[10]]$mm)
    vv<-cbind(all.data[[1]]$vv, all.data[[2]]$vv, all.data[[3]]$vv, all.data[[4]]$vv, all.data[[5]]$vv, all.data[[6]]$vv, all.data[[7]]$vv, all.data[[8]]$vv, all.data[[9]]$vv, all.data[[10]]$vv)
    nn<-cbind(all.data[[1]]$nn, all.data[[2]]$nn, all.data[[3]]$nn, all.data[[4]]$nn, all.data[[5]]$nn, all.data[[6]]$nn, all.data[[7]]$nn, all.data[[8]]$nn, all.data[[9]]$nn, all.data[[10]]$nn)
    tt<-cbind(all.data[[1]]$tt, all.data[[2]]$tt, all.data[[3]]$tt, all.data[[4]]$tt, all.data[[5]]$tt, all.data[[6]]$tt, all.data[[7]]$tt, all.data[[8]]$tt, all.data[[9]]$tt, all.data[[10]]$tt)
  }

  y <- list(xx = xx, vv = vv, nn = nn, tt = tt)

  return(y)
}
