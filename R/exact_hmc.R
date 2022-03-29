#' Title
#'
#' @param 
#' @return
#' @export
#'
#' @examples
run_bouncy_sampler = function(n,
                              initial_position,
                              constraint_direc,
                              constraint_bound,
                              cholesky,
                              mu,
                              precision = TRUE,
                              total_time = pi / 2,
                              seed = 1) {
  set.seed(seed)
  results = matrix(nrow = ncol(constraint_direc), ncol = n)
  whitened_constraints = WhitenConstraints(constraint_direc,
                                           constraint_bound,
                                           cholesky,
                                           mu,
                                           precision)
  sample = initial_position
  for (i in 1:n) {
    initial_momentum = rnorm(ncol(constraint_direc))
    sample = GenerateSample(
      sample,
      initial_momentum,
      whitened_constraints$direc,
      whitened_constraints$direc_rownorm_sq,
      whitened_constraints$bound,
      cholesky,
      mu,
      total_time,
      precision
    )
    results[, i] = sample
  }
  return(results)
}