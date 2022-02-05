#' Run one iteration of the standard HMC on a constrained Gaussian target
#'
#' As opposed to zigzag HMC that relies on Laplace momentum, the standard HMC 
#' deploys a Gaussian momentum. When the target is Gaussian, the resulting 
#' dynamics consists of independent harmonic oscillations along principal components.
#'
#' @param cholFactor
#'   Cholesky factor of the precision or covariance as specified by `paramet`
#' @param paramet
#'   string 'prec' or 'cov', specifying whether the Cholesky factor corresponds
#'   to that of the precision or covariance matrix
harmonicHmc <- function(
  mean, cholFactor, constraitDirec, constraintBound, t, paramet = 'prec'
) {
  
}

#' Sample from a constrained Guassian using harmonic HMC
#' 
#' @param n number of samples
#' @param mean 
#'   d-dimensional vector of the pre-truncation mean of the Gaussian distribution
#' @param prec 
#'   d-by-d precision matrix of the Gaussian distribution
#' @param cov 
#'   d-by-d covariance matrix of the Gaussian distribution; exactly one of 
#'   `prec` or `cov` must be specified
#' @param constraintDirec
#'   k-by-d matrix where k is the number of linear constraints
#' @param constraintBound
#'   k-dimensional vector. `constraintDirec` and `constraintBound` together 
#'   specify the constraint `constraintDirec %*% param >= constraintBound`
#' 
harmonicHmcSampler <- function(
  n, mean, cov = NULL, prec = NULL, constraitDirec, constraintBound, t, 
  burnin = 0, p0 = NULL, random_seed = 666
) {
  
}

#' Sample from a constrained Guassian using harmonic HMC
#' 
#' @param constraints
#'   vector of +1 or -1 indicating orthants, with 0 indicating no constraint 
#'   along the coordinate (?)
harmonicHmcSamplerWrapper <- function(
  n, mean, cov = NULL, prec = NULL, constraits, t, 
  burnin = 0, p0 = NULL, random_seed = 666
) {
  
}