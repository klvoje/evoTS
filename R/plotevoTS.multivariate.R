#'@title Plots multivariate evolutionary sequence (time-series) data set
#'
#'@description Function to plot multivariate evolutionary sequence (time-series), showing trait means over time.
#'
#' @param yy a multivariate evoTS object
#'
#' @param nse the number of standard errors represented by the error bars on the plot; default is 1
#'
#' @param col vector indicating colors
#'
#' @param lty line type
#'
#' @param lwd line width
#'
#' @param pch plotting symbols
#'
#' @param x.label label on x axis
#'
#' @param y.label label on y axis
#'
#' @param y_min minimum value of y axis
#'
#' @param y_max maximum value of y axis
#'
#' @param cex.axis Specify the size of the tick label numbers/text
#'
#' @param cex.lab specify the size of the axis label
#'
#' @param cex.main specify the size of the title text
#'
#' @param axes logical, whether to plot axes or not
#'
#'@return The results are plotted.
#'
#'@author Kjetil Lysne Voje
#'
#'@export
#
#'@examples
#'## Generate two evolutionary sequences (time-series)
#'x1 <- paleoTS::sim.Stasis(60, vp=1)
#'x2 <- paleoTS::sim.Stasis(60, vp=1)
#'
#'## Make a multivariate data set
#'x1_x2<-make.multivar.evoTS(x1, x2)
#'
#'## Plot the data
#'plotevoTS.multivariate(x1_x2, y_min=-1, y_max=1)


plotevoTS.multivariate<-function(yy, nse = 1, col = NULL, lty = NULL, lwd = NULL, pch = NULL, x.label=NULL, y.label=NULL, y_min=NULL, y_max=NULL, cex.axis=NULL, cex.lab=NULL, cex.main=NULL, axes=NULL){

  m <-ncol(yy$xx) # number of traits
  if (is.null(col)) {
    col<-c("cornflowerblue", "orange", "red2", "forestgreen", "black", "darkblue", "tan4", "indianred", "mediumorchid", "firebrick4")
  }

  if (is.null(lty)) lty=1
  if (is.null(lwd)) lwd=1
  if (is.null(pch)) pch=1
  if (is.null(x.label)) x.label="Time"
  if (is.null(y.label)) y.label="Trait mean"
  if (is.null(y_min))   y_min<-min(yy$xx[,])
  if (is.null(y_max))   y_max<-max(yy$xx[,])
  if (is.null(cex.axis))   cex.axis<-1
  if (is.null(cex.lab))   cex.lab<-1
  if (is.null(cex.main))   cex.main<-1
  if (is.null(axes))   axes<-TRUE

  SE.1<-sqrt(yy$vv[,1]/yy$nn[,1])*nse
  SE.2<-sqrt(yy$vv[,2]/yy$nn[,2])*nse

  plot(yy$tt[,1],yy$xx[,1],col=col[1], ylim=c(y_min,y_max),  pch=pch, xlab=x.label, ylab=y.label, cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main, axes = axes)
  lines(yy$tt[,1], yy$xx[,1], col = col[1])
  graphics::segments(x0= yy$tt[,1], y0=yy$xx[,1]+SE.1, x1=yy$tt[,1], y1=yy$xx[,1]-SE.1, col = col[1])
  points(yy$tt[,2],yy$xx[,2], col=col[2], pch=pch)
  lines(yy$tt[,2],yy$xx[,2], col=col[2])
  graphics::segments(x0= yy$tt[,2], y0=yy$xx[,2]+SE.2, x1=yy$tt[,2], y1=yy$xx[,2]-SE.2, col = col[2])


  if (m>2){
    for (i in 3:m){
      points(yy$tt[,i],yy$xx[,i], col=col[i])
      lines(yy$tt[,i],yy$xx[,i], col=col[i])
      graphics::segments(x0= yy$tt[,i], y0=yy$xx[,i]+(sqrt(yy$vv[,i]/yy$nn[,i])*nse), x1=yy$tt[,i], y1=yy$xx[,i]-(sqrt(yy$vv[,i]/yy$nn[,i])*nse), col = col[i])

    }
  }

}
