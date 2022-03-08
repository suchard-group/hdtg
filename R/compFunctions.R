#' Title
#'
#' @param X 
#' @param n 
#' @param p 
#' @param sigma2 
#'
#' @return
#' @export
#'
#' @examples
getPrecCov <- function(X, n, p, sigma2 = 1){
  
  A = 1 / sigma2 * diag(p) + t(X) %*% X
  B = -t(X)
  C = -X
  D = diag(n)
  precMat = rbind(cbind(A, B), cbind(C, D))
  
  A_cov = sigma2 * diag(p)
  B_cov = sigma2 * t(X)
  C_cov = sigma2 * X
  D_cov = diag(n) + sigma2 * X %*% t(X)
  covMat = rbind(cbind(A_cov, B_cov), cbind(C_cov, D_cov))
  assertthat::assert_that(!any(is.na(covMat)))
  return(list(precMat = precMat, covMat = covMat))
}