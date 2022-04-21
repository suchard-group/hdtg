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
#' @param total_time amount of time the particle bounces for each sample
#' @param seed random seed
#' @param diagnostic_mode boolean for whether to return the bounce distances for
#' each sample
#' @return List of
#' "samples": d x n matrix of samples
#' "bounces_distances": list of bounces for each sample
#' @export
#'
#' @examples
#' set.seed(1)
#' d = 100
#' A = matrix(runif(d^2)*2-1, ncol=d)
#' Sigma = t(A) %*% A
#' R = cholesky(solve(Sigma))
#' mu = rep(0,d)
#' constraint_direc = diag(d)
#' constraint_bound = rep(0,d)
#' initial = rep(1, d)
#' results = runBouncySampler(
#' 10000,
#' initial,
#' constraint_direc,
#' constraint_bound,
#' R,
#' mu,
#' )

runBouncySampler = function(n,
                            initial_position,
                            constraint_direc,
                            constraint_bound,
                            cholesky_factor,
                            unconstrained_mean,
                            prec_parametrized = TRUE,
                            total_time = pi / 2,
                            seed = 1,
                            diagnostic_mode = FALSE) {
  set.seed(seed)
  samples = matrix(nrow = ncol(constraint_direc), ncol = n)
  bounce_distances = vector(mode = "list",
                            length = ifelse(diagnostic_mode, n, 0))
  whitened_constraints = applyWhitenTransform(
    constraint_direc,
    constraint_bound,
    cholesky_factor,
    unconstrained_mean,
    prec_parametrized
  )
  sample = initial_position
  for (i in 1:n) {
    initial_momentum = rnorm(ncol(constraint_direc))
    results = generateSample(
      sample,
      initial_momentum,
      whitened_constraints$direc,
      whitened_constraints$direc_rownorm_sq,
      whitened_constraints$bound,
      cholesky_factor,
      unconstrained_mean,
      total_time,
      prec_parametrized,
      diagnostic_mode
    )
    sample = results$sample
    samples[, i] = sample
    if (diagnostic_mode) {
      bounce_distances[[i]] = results$bounce_distances
    }
  }
  return(list("samples" = samples, "bounce_distances" = bounce_distances))
}
