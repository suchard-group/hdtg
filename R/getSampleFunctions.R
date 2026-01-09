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

#' Draw one Markovian zigzag sample
#'
#' Simulate the Markovian zigzag dynamics for a given position over a specified travel time.
#'
#' @param position a d-dimensional position vector.
#' @param velocity optional d-dimensional velocity vector. If NULL, it will be generated within the function.
#' @param engine an object representing the Markovian zigzag engine, typically containing settings and state required for the simulation.
#' @param travelTime the duration for which the dynamics are simulated.
#'
#' @return A list containing the position and velocity after simulating the dynamics.
#' @export
getMarkovianZigzagSample <- function(position,
                                     velocity = NULL,
                                     engine,
                                     travelTime) {
  if (is.null(velocity)) {
    velocity <- 2 * stats::rbinom(length(position), 1, .5) - 1
  }
  
  res <- .oneIrreversibleIteration(
    sexp = engine$engine,
    position = position,
    velocity = velocity,
    time = travelTime
  )
  return(list(position=res$position, velocity=res$velocity))
}