#' Sample from a constrained Gaussian distribution
#'
#' Generate samples from a constrained Gaussian distribution with the Hamiltonian zigzag sampler.
#'
#' @param n number of samples
#' @param mean a d-dimensional mean vector
#' @param prec a d-by-d precision matrix of the Gaussian distribution
#' @param constraits a d-dimensional vector where a positive (negative) number means >0 (<0) truncation.
#' @param p0 a d-dimensional vector of the initial value. It must satisfy all constraints. If not specified a random initial value will be used
#' @param burnin the number of burn-in iterations
#' @param t time length to simulate the Markov process
#'
#' @return An n-by-d matrix where each row is a multivariate sample
#' @export
#'
#' @examples rcmg(1,1,1,1,1)
rcmg <- function(n, mean, prec, constraits, t, burnin, p0 = NULL) {
  debug_flg = T
  stopifnot("n > burnin must be integers!" = (n %% 1 == 0 && burnin %% 1 == 0 && n > burnin))
  stopifnot("mean and prec must be numeric" = (is.numeric(mean) && is.numeric(prec)))
  # TODO add other checks for arguments. all dimensions must match.
  
  ndim = length(mean)
  get_prec_product <- function (x) {
    if (length(x) == 1){
      return(prec[, x])
    } else {
      return(drop(prec %*% x))
    }
  }
  
  samples <- array(0, c(ndim, n))
  set.seed(666)

  #
  for (i in 1:n) {
    momentum <-
      (2 * (runif(ndim) > .5) - 1) * rexp(ndim, rate = 1)
    t_jittered <- t
    p0 <- hzz(get_prec_product, mean, p0, constraits, momentum, t_jittered)
    samples[, i] <- p0
    if (debug_flg) {
      cat('iteration', i, 'done \n')
    }
  }
  return(samples)
}
