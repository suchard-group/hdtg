#' Sample from a constrained Gaussian distribution
#'
#' Generate samples from a constrained Gaussian distribution with the Hamiltonian zigzag sampler.
#'
#' @param n number of samples
#' @param mean a d-dimensional mean vector
#' @param prec a d-by-d precision matrix of the Gaussian distribution
#' @param constraits a (list?) of constraints
#' @param p0 a d-dimensional vector of the initial value. It must satisfy all constraints. If not specified a random initial value will be used
#' @param burnin the number of burn-in iterations
#'
#' @return An n-by-d matrix where each row is a multivariate sample
#' @export
#'
#' @examples rcmg(1,1,1,1,1)
rcmg <- function(n, mean, prec, constraits, p0, burnin) {
  return(0)
}
