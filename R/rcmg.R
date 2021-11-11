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
#' @param t time length to simulate the Markov process
#'
#' @return An n-by-d matrix where each row is a multivariate sample
#' @export
#'
#' @examples rcmg(1,1,1,1,1)
rcmg <- function(n, mean, prec, constraits, t, burnin, p0 = NULL) {
  
  stopifnot("n > burnin must be integers!" = (n %% 1 == 0 && burnin %% 1 == 0 && n > burnin))
  stopifnot("mean and prec must be numeric" = (is.numeric(mean) && is.numeric(prec)))
  # TODO add other checks for arguments. all dimensions must match.
  
  ndim = length(mean)
  get_prec_product <- function (x) {
      drop(prec %*% x)
  }
  
  samples <- array(0, c(ndim, n))
  for (i in 1:n) {
    t_jittered <- t + .1 * runif(1, -t, t)
    
    x <- hzz(get_prec_product, mean, position, momentum, t_jittered)

    samples[, i] <- x
  }
  
 
  
  return(samples)
}
