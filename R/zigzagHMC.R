#' Sample from a truncated Gaussian distribution 
#'
#' Generate MCMC samples from a d-dimensional truncated Gaussian distribution with element-wise truncations using the Zigzag Hamiltonian Monte Carlo sampler (Zigzag-HMC).
#'
#' @param nSample number of samples after burn-in.
#' @param burnin number of burn-in samples (default = 0).
#' @param mean a d-dimensional mean vector.
#' @param prec a d-by-d precision matrix of the Gaussian distribution. 
#' @param lowerBounds a d-dimensional vector specifying the lower bounds. `-Inf` is accepted.  
#' @param upperBounds a d-dimensional vector specifying the upper bounds. `Inf` is accepted. 
#' @param nutsFlg logical. If `TRUE` the No-U-Turn sampler will be used (Zigzag-NUTS).
#' @param precondition logical. If `TRUE`, the precision matrix will be preconditioned so that its diagonals (i.e. conditional variances) are all 1.
#' @param init a d-dimensional vector of the initial value. `init` must satisfy all constraints. If `init = NULL`, a random initial value will be used.
#' @param stepSize step size for Zigzag-HMC or Zigzag-NUTS (if `nutsFlg = TRUE`). Default value is the empirically optimal choice: sqrt(2)(lambda)^(-1/2) for Zigzag-HMC and 0.1(lambda)^(-1/2) for Zigzag-NUTS, where lambda is the minimal eigenvalue of the precision matrix.   
#' @param seed random seed (default = 1).
#' @param numThreads number of threads for parallel execution (default = 1). Set to 0 for automatic detection of available cores.
#' @param diagnosticMode logical. `TRUE` for also returning diagnostic information such as the stepsize used. 
#' @return 
#' When `diagnosticMode = FALSE` (default), returns an `nSample`-by-`d` matrix of samples.
#' 
#' When `diagnosticMode = TRUE`, returns a list with elements:
#' \item{samples}{`nSample`-by-`d` matrix of samples}
#' \item{stepsize}{The step size used for sampling}
#' @export
#' @examples
#' set.seed(1)
#' d <- 10
#' A <- matrix(runif(d^2)*2-1, ncol=d)
#' precMat <- t(A) %*% A
#' initial <- rep(1, d)
#' results <- zigzagHMC(nSample = 1000, burnin = 1000, mean = rep(0, d), prec = precMat,
#' lowerBounds = rep(0, d), upperBounds = rep(Inf, d))
#'
#' @references
#'
#' Nishimura, A., Zhang, Z., and Suchard, M. A. (2024). Zigzag path connects 
#' two Monte Carlo samplers: Hamiltonian counterpart to a piecewise deterministic 
#' Markov process. Journal of the American Statistical Association, 1-13.
#'
#' Nishimura, A., Dunson, D. B., and Lu, J. (2020). Discontinuous Hamiltonian 
#' Monte Carlo for discrete parameters and discontinuous likelihoods. 
#' Biometrika, 107(2): 365-380.
#' 
#' @seealso [getZigzagSample()], [createEngine()], [createNutsEngine()], [setMean()], [setPrecision()] 
zigzagHMC <- function(nSample,
                      burnin = 0,
                      mean,
                      prec,
                      lowerBounds,
                      upperBounds,
                      init = NULL,
                      stepSize = NULL,
                      nutsFlg = FALSE,
                      precondition = FALSE,
                      seed = 1,
                      numThreads = 1,
                      diagnosticMode = FALSE) {
  
  validateInput(mean, prec, lowerBounds, upperBounds, init)
  if (is.null(init)) {
    init <- getInitialPosition(mean, lowerBounds, upperBounds)
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  cpp_seed <- sample.int(.Machine$integer.max, size = 1)
  
  if (precondition) {
    precondScaleFactor <- sqrt(diag(prec))
    init <- precondScaleFactor * init
    mean <- precondScaleFactor * mean
    prec <- stats::cov2cor(prec)
  }
  
  ndim <- length(mean)
  samples <- array(0, c(nSample, ndim))
  
  if (nutsFlg) {
    if (is.null(stepSize)) {
      stepSize <- 0.1 / sqrt(computeExtremeEigenval(prec))
    }
    engine <- createNutsEngine(
      dimension = ndim,
      lowerBounds = lowerBounds,
      upperBounds = upperBounds,
      flags = 128,
      numThreads = numThreads,
      seed = cpp_seed,
      stepSize = stepSize,
      mean = mean,
      precision = prec
    )
    
  } else {
    if (is.null(stepSize)) {
      stepSize <- sqrt(2) / sqrt(computeExtremeEigenval(prec))
    }
    engine <- createEngine(
      dimension = ndim,
      lowerBounds = lowerBounds,
      upperBounds = upperBounds,
      flags = 128,
      numThreads = numThreads,
      seed = cpp_seed,
      mean = mean,
      precision = prec
    )
  }
  
  position <- init
  for (i in 1:(nSample + burnin)) {
    position <- getZigzagSample(
      position = position,
      momentum = NULL,
      nutsFlg = nutsFlg,
      engine = engine,
      stepSize = stepSize
    )
    if (i > burnin) {
      if (precondition) {
        samples[i - burnin, ] <- position / precondScaleFactor
      } else {
        samples[i - burnin, ] <- position
      }
    }
  }
  if (diagnosticMode) {
    return(list("samples" = samples, "stepsize" = stepSize))
  } else {
    return(samples)
  }
}

#' Draw one MTN sample with Zigzag-HMC or Zigzag-NUTS
#'
#' Simulate the Zigzag-HMC or Zigzag-NUTS dynamics on a given MTN.
#'
#' @param position a d-dimensional initial position vector.
#' @param momentum a d-dimensional initial momentum vector.
#' @param nutsFlg logical. If `TRUE` the No-U-Turn sampler will be used (Zigzag-NUTS).
#' @param engine list. Its `engine` element is a pointer to the Zigzag-HMC engine
#' (or Zigzag-NUTS engine) C++ object that implements fast computations for
#' Zigzag-HMC (or Zigzag-NUTS).
#' @param stepSize step size for Zigzag-HMC. If `nutsFlg = TRUE`, `engine` contains
#' the base step size for Zigzag-NUTS).
#'
#' @return one MCMC sample from the target MTN.
#' @export
#' @note `getZigzagSample` is particularly efficient when the target MTN has a random
#' mean and covariance/precision where one can reuse the Zigzag-HMC engine object while
#' updating the mean and covariance. The following example demonstrates such a use.

#' @examples 
#' set.seed(1)
#' n <- 1000
#' d <- 10
#' samples <- array(0, c(n, d))
#' 
#' # initialize MTN mean and precision
#' m <- rnorm(d, 0, 1)
#' prec <- rWishart(n = 1, df = d, Sigma = diag(d))[, , 1]
#' # call createEngine once
#'engine <- createEngine(dimension = d, lowerBounds = rep(0, d),
#'  upperBounds = rep(Inf, d), seed = 1, mean = m, precision = prec)
#'
#' HZZtime <- sqrt(2) / sqrt(min(mgcv::slanczos(
#'  A = prec, k = 1,
#'  kl = 1
#' )[['values']]))
#'
#' currentSample <- rep(0.1, d)
#' for (i in 1:n) {
#'   m <- rnorm(d, 0, 1)
#'   prec <- rWishart(n = 1, df = d, Sigma = diag(d))[,,1]
#'   setMean(engine = engine, mean = m)
#'   setPrecision(engine = engine, precision = prec)
#'   currentSample <- getZigzagSample(position = currentSample,
#'                                    nutsFlg = FALSE,
#'                                    engine = engine,
#'                                    stepSize = HZZtime)
#'   samples[i,] <- currentSample
#' }
#' @seealso [zigzagHMC()], [drawLaplaceMomentum()]
getZigzagSample <- function(position,
                            momentum = NULL,
                            nutsFlg,
                            engine,
                            stepSize = NULL) {
  if (is.null(momentum)) {
    momentum <- drawLaplaceMomentum(length(position))
  }
  
  if (nutsFlg) {
    res <- .oneNutsIteration(sexp = engine$engine,
                             position = position,
                             momentum = momentum)
    
  } else {
    res <- .oneIteration(
      sexp = engine$engine,
      position = position,
      momentum = momentum,
      time = stepSize
    )
  }
  return(res$position)
}