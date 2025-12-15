#' @useDynLib hdtg, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom stats rnorm runif
#' @importFrom Rdpack reprompt
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("hdtg", libpath)
}
