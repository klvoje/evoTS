#' Evolutionary sequence (time-series) of phenotypic change in diameter in the lineage Stephanodiscus yellowstonensis
#'
#' Phenotypic data (diameter) from a centric diatom lineage Stephanodiscus yellowstonensis. The time series spans about 14 000 years.
#' The data set contains data on valve diameter (measured in micrometres).
#' The data consists of an object of class paleoTS (diameter_S.yellowstonensis). Objects of class paleoTS can be analyzed in evoTS.
#' The object (trait data set) contains a vector of sample means (mm), sample variances (vv), sample sizes (nn) and sample ages (tt).
#' The oldest sample is listed first. The data spans an interval of 13728 years.
#'
#' @docType data
#'
#' @usage data(diameter_S.yellowstonensis)
#'
#' @format An object of class \code{"paleoTS"}.
#'
#' @keywords datasets
#'
#' @references Theriot et al. 2006. Late Quaternary rapid morphological evolution of an endemic diatom in Yellowstone Lake, Wyoming. Paleobiology 32:38-54 
#'
#' @examples
#' ln.diameter<-paleoTS::ln.paleoTS(diameter_S.yellowstonensis)
#' ln.diameter$tt<-ln.diameter$tt/(max(ln.diameter$tt))
#' opt.joint.decel(ln.diameter)
"diameter_S.yellowstonensis"
