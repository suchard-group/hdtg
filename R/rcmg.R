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
#' @param cpp_flg 
#'
#' @return An n-by-d matrix where each row is a multivariate sample
#' @export
#'
#' @examples rcmg(1,1,1,1,1)
rcmg <- function(n, mean, cov = NULL, prec = NULL, constraits, t, burnin, p0 = NULL, cpp_flg = FALSE, nuts_flg = FALSE, random_seed = 666, debug_flg = F) {
  stopifnot("n > burnin must be integers!" = (n %% 1 == 0 && burnin %% 1 == 0 && n > burnin))
  stopifnot("mean and prec must be numeric" = (is.numeric(mean) && is.numeric(prec)))
  stopifnot("must provide either covariance or precision" = (!is.null(cov) || !is.null(prec)))
  if (is.null(cov)){
    cov <- solve(prec)
  }
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
  set.seed(random_seed)

  if (cpp_flg) {
    engine <- createEngine(dimension = ndim, mask = rep(1, ndim), observed = rep(1, ndim), parameterSign = constraits, flags = 128, info = 1, seed = 1, mean = mean, precision = prec)
  } else {
    engine <- NULL
  }

  for (i in 1:n) {
    momentum <-
      (2 * (runif(ndim) > .5) - 1) * rexp(ndim, rate = 1)
    t_jittered <- t
    
    p0 <- hzz(get_prec_product = get_prec_product, mean = mean, position = p0, constraits = constraits, momentum = momentum, t = t_jittered, cpp_flg = cpp_flg, nuts_flg = nuts_flg, engine = engine)
    
    samples[, i] <- p0
    if (debug_flg) {
      cat('iteration', i, 'done \n')
    }
  }
  return(samples)
}

