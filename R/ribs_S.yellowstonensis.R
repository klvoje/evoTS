#' Evolutionary sequence (time-series) of phenotypic change in the number of ribs in the lineage Stephanodiscus yellowstonensis
#'
#' Phenotypic data (ribs) from a centric diatom lineage Stephanodiscus yellowstonensis. The time series spans about 14 000 years.
#' The data set contains data on the number of costae (ribs) per valve.
#' The data consists of an objects of class paleoTS (ribs_S.yellowstonensis). Objects of class paleoTS can be analyzed in evoTS.
#' The object (trait data set) contains a vector of sample means (mm), sample variances (vv), sample sizes (nn) and sample ages (tt).
#' The oldest sample is listed first. The data spans an interval of 13728 years.
#'
#' @docType data
#'
#' @usage data(ribs_S.yellowstonensis)
#'
#' @format An object of class \code{"paleoTS"}.
#'
#' @keywords datasets
#'
#' @references Theriot et al. 2006. Late Quaternary rapid morphological evolution of an endemic diatom in Yellowstone Lake, Wyoming. Paleobiology 32:38-54 
#'
#' @examples
#' ln.ribs<-paleoTS::ln.paleoTS(ribs_S.yellowstonensis)
#' ln.ribs$tt<-ln.ribs$tt/(max(ln.ribs$tt))
#' opt.joint.decel(ln.ribs)
"ribs_S.yellowstonensis"
