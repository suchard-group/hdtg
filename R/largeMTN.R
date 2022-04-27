#' @useDynLib largeMTN, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom RcppXsimd supportsSSE supportsAVX
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("largeMTN", libpath)
}