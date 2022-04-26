#' @useDynLib largeMTN, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("largeMTN", libpath)
}