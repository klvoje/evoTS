#'@title Plot a paleoTS object 
#'
#'@description Plot a paleoTS object (slightly modified version of the same function in paleoTS)
#'
#' @param x a paleoTS object
#'
#' @param nse the number of standard errors represented by the error bars on the plot; defaults to 1
#'
#' @param pool logical indicating if variances should be pooled across samples for the purposes of displaying error bars; defaults to FALSE
#'
#' @param add logical, if TRUE, adds to existing plot
#'
#' @param modelFit optional model fit from fitting functions
#'
#' @param pch plotting symbol, defaults to 19
#'
#' @param lwd line width, defaults to 1.5
#'
#' @param ylim optional, y-limits of the plot
#' 
#' @param xlab a title for the x axis
#' 
#' @param ylab a title for the y axis
#'
#' @param ... other arguments passed to plotting functions
#'
#'@return The results are plotted.
#'
#'@author Kjetil Lysne Voje
#'
#'@export
#
plotevoTS <-function (x, nse = 1, pool = FALSE, add = FALSE, modelFit = NULL, 
                        pch = 19, lwd = 1.5, ylim=NULL, xlab=NULL, ylab=NULL, ...) 
{
  if (pool) 
    x <- paleoTS::pool.var(x, ret.paleoTS = TRUE)
  se <- sqrt(x$vv/x$nn)
  lci <- x$mm - (nse * se)
  uci <- x$mm + (nse * se)
  xx <- x
  if (!is.null(x$start.age)) {
    if(x$timeDir=="decreasing" & max(x$tt)==1) {x$tt <- x$tt; xl <- (range(x$tt))}
    if(x$timeDir=="decreasing" & max(x$tt)!=1) { x$tt <- x$start.age - x$tt; xl <- rev(range(x$tt))}
    if(x$timeDir=="increasing")	{x$tt<- x$tt + x$start.age; xl<- range(x$tt)}
  }
  
  else xl <- range(x$tt)
  
  if(is.null(xlab)) xlab<- "Time"
  if(is.null(ylab)) ylab<- "Trait Mean"
  
  
  if (!is.null(modelFit)) {
    mlab <- paste(modelFit$modelName, "expectation [95% prob. interval]")
    mc <- modelCurves(xx, w = modelFit)
    if (is.na(mc$ee[1])) 
      modelFit <- NULL
  }
  if (is.null(modelFit)) 
    yl <- c(uci, lci)
  else yl <- c(uci, lci, mc$ll, mc$uu)
  if(is.null(ylim)) ylim<- range(yl)
  if (!add) 
    plot(range(x$tt), ylim=ylim, typ = "n", pch = 19, 
         xlab = xlab, ylab = ylab, xlim = xl, ...)
  if (!is.null(modelFit)) {
    if (!is.null(x$start.age)) 
      mc$tt <- x$start.age - mc$tt
    polygon(c(mc$tt, rev(mc$tt)), c(mc$uu, rev(mc$ll)), col = "wheat2", 
            border = "white")
    lines(mc$tt, mc$ee, col = "tan", lwd = 2)
  }
  lines(x$tt, x$mm, lwd = lwd, ...)
  segments(x$tt, lci, x$tt, uci, lty = 1, lwd = lwd, ...)
  points(x$tt, x$mm, pch = pch, cex = 1.2, ...)
  mtext(x$label, cex = 0.7, col = "grey", font = 3)
  if (!is.null(modelFit)) 
    mtext(mlab, side = 4, cex = 0.8, col = "tan", font = 2)
}