#' @useDynLib ZigZag, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("mypackage", libpath)
}