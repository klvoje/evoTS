#' @title Fit univariate punctuation models to an evolutionary sequence (time-series).
#'
#' @description Wrapper function to find maximum likelihood solutions for univariate univariate punctuation models to an evolutionary sequence (time-series).
#'
#' @param y an univariate paleoTS object.
#' 
#' @param minb the minimum number of samples within a segment to consider
#' 
#' @param fit.2.punc logical, indicating whether to fit models with 2 mode shifts
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
#'## Fit univariate punctuation models (the example may take > 5 seconds to run)
#'##fit.PE.models(x)
#'

fit.PE.models<-function (y, minb=7, fit.2.punc=FALSE)
{
  y$start.age<-NULL
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
  
  if (fit.2.punc==FALSE){
    m1 <- paleoTS::fitGpunc(y, method="Joint", minb = minb, oshare=FALSE)
    m2 <- fit.Bokma(y, minb = minb)
  }

  if (fit.2.punc==TRUE){
    m1 <- paleoTS::fitGpunc(y, method="Joint", minb = minb, ng = 3, oshare=FALSE)
    m2 <- fit.2.Bokma(y, minb = minb)
  }
  
  mc <- paleoTS::compareModels(m1, m2, silent = FALSE)
  invisible(mc)
  
}
