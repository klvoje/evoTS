#' @title Fit all univariate models to an evolutionary sequence (time-series).
#'
#' @description Wrapper function to find maximum likelihood solutions for all univariate models (excluding models with mode shifts) to an evolutionary sequence (time-series).
#'
#' @param y an univariate paleoTS object.
#' 
#' @param pool indicating whether to pool variances across samples
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
#'## ##Generate a paleoTS object.
#'x <- paleoTS::sim.GRW(30)
#'
#'## Fit univariate models to the data.
#'fit.all.univariate(x, pool = TRUE)
#'

fit.all.univariate<-function (y, pool = TRUE)
{
  y$start.age<-NULL

  if (pool) 
    y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  
    m1 <- paleoTS::opt.joint.GRW(y, pool = pool)
    m2 <- paleoTS::opt.joint.URW(y, pool = pool)
    m3 <- paleoTS::opt.joint.Stasis(y, pool = pool)
    m4 <- opt.joint.StrictStasis(y, pool = pool)
    m5 <- opt.joint.decel(y, pool = pool)
    m6 <- opt.joint.accel(y, pool = pool)
    m7 <- paleoTS::opt.joint.OU(y, pool = pool)
    m8 <- opt.joint.OUBM(y, opt.anc = TRUE, pool = pool)
    m9 <- opt.joint.OUBM(y, opt.anc = FALSE, pool = pool)

  mc <- paleoTS::compareModels(m1, m2, m3, m4, m5, m6, m7, m8, m9, silent = FALSE)
  invisible(mc)
}
