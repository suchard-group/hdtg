#' Run Bouncy HMC Sampler
#' 
#' Sample from a truncated Gaussian distribution with constraints Fx+g >= 0.
#'
#' @param n number of samples
#' @param initial_position starting value for parameters
#' @param constraint_direc F matrix (k-by-d matrix where k is the number of 
#' linear constraints)
#' @param constraint_bound g vector (k dimensional)#' 
#' @param cholesky_factor upper triangular matrix R from cholesky decomposition of 
#' precision or covariance matrix into R^TR
#' @param unconstrained_mean mean of unconstrained Gaussian
#' @param prec_parametrized boolean for whether parametrization is by precision (TRUE) 
#' or covariance matrix (FALSE)
#' @param total_time amount of time the particl bounces for each sample
#' @param seed random seed
#' @return d x n matrix of samples 
#' @export
#'
#' @examples
#' set.seed(1)
#' d = 100
#' A = matrix(runif(d^2)*2-1, ncol=d)
#' Sigma = t(A) %*% A
#' R = Cholesky(solve(Sigma))
#' mu = rep(0,d)
#' constraint_direc = diag(d)
#' constraint_bound = rep(0,d)
#' results = run_bouncy_sampler(
#' 10000,
#' rep(1, d),
#' constraint_direc,
#' constraint_bound,
#' R,
#' mu,
#' )

run_bouncy_sampler = function(n,
                              initial_position,
                              constraint_direc,
                              constraint_bound,
                              cholesky_factor,
                              unconstrained_mean,
                              prec_parametrized = TRUE,
                              total_time = pi / 2,
                              seed = 1) {
  set.seed(seed)
  results = matrix(nrow = ncol(constraint_direc), ncol = n)
  whitened_constraints = ApplyWhitenTransform(constraint_direc,
                                              constraint_bound,
                                              cholesky_factor,
                                              unconstrained_mean,
                                              prec_parametrized)
  sample = initial_position
  for (i in 1:n) {
    initial_momentum = rnorm(ncol(constraint_direc))
    sample = GenerateSample(
      sample,
      initial_momentum,
      whitened_constraints$direc,
      whitened_constraints$direc_rownorm_sq,
      whitened_constraints$bound,
      cholesky_factor,
      unconstrained_mean,
      total_time,
      prec_parametrized
    )
    results[, i] = sample
  }
  return(results)
}
