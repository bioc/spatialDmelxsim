#' Spatial allelic expression counts for fly cross embryo
#'
#' Allelic expression counts of spatial slices of a fly embryo from
#' Combs & Fraser (2018), a D melanogaster x D simulans reciprocal cross
#'
#' @references
#'
#' Combs PA, Fraser HB (2018) Spatially varying cis-regulatory divergence
#' in Drosophila embryos elucidates cis-regulatory logic.
#' PLOS Genetics 14(11): e1007631.
#'
#' @importFrom utils read.csv
#' @docType package
#' @name spatialDmelxsim
NULL

.onLoad <- function(libname, pkgname) {
    fl <- system.file("extdata", "metadata.csv", package=pkgname)
    titles <- read.csv(fl, stringsAsFactors=FALSE)$Title
    ExperimentHub::createHubAccessors(pkgname, titles)
}
