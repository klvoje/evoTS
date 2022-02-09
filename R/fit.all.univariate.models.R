#' @title Fit all univariate models to an evolutionary sequence (time-series).
#'
#' @description Wrapper function to find maximum likelihood solutions for all univariate models to an evolutionary sequence (time-series).
#'
#' @param y an univariate paleoTS object.
#'
#'@return The function returns a list of all investigated models and their highest log-likelihood (and their corresponding AICc and AICc weight).
#'
#'@author Kjetil Lysne Voje
#'
#'@references Hunt, G. 2006. Fitting and comparing models of phyletic evolution: random walks and beyond. \emph{Paleobiology} 32:578–601
#'@references Hunt, G., Bell, M. A. & Travis, M. P. Evolution towards a new adaptive optimum: Phenotypic evolution in a fossil stickleback lineage. \emph{Evolution} 62:700–710 (2008)
#'
#'@export
#'
#'@examples
#'## Generate an univariate evolutionary sequence data set
#'x <- paleoTS::sim.GRW(30)
#'
#'## Fit a the model to the data by defining shift points.
#'fit.all.univariate.models(x)
#'

fit.all.univariate.models<-function (y)
{
  args <- list()
  check.var <- TRUE
  if (length(args) > 0)
    if (args$pool == FALSE)
      check.var <- FALSE
  if (check.var) {
    tv <- paleoTS::test.var.het(y)
    pv <- round(tv$p.value, 0)
    wm <- paste("Sample variances not equal (P = ", pv, "); consider using argument pool=FALSE",
                collapse = "")
    if (pv <= 0.05)
      warning(wm)
  }
    m1 <- paleoTS::opt.joint.GRW(y)
    m2 <- paleoTS::opt.joint.URW(y)
    m3 <- paleoTS::opt.joint.Stasis(y)
    m4 <- opt.joint.StrictStasis(y)
    m5 <- opt.joint.decel(y)
    m6 <- opt.joint.accel(y)
    m7 <- paleoTS::opt.joint.OU(y)
    m8 <- opt.joint.OUBM(y, opt.anc = TRUE)
    m9 <- opt.joint.OUBM(y, opt.anc = FALSE)
    m10 <- paleoTS::fitGpunc(y, method="Joint")
    m11 <- fit.Bokma(y)

  mc <- paleoTS::compareModels(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, silent = FALSE)
  invisible(mc)
}
