#' Markovian Zigzag Sampler
#'
#' Sample from a truncated multivariate normal distribution using the
#' Markovian Zigzag process, a continuous-time, non-reversible Markov chain
#' Monte Carlo method based on piecewise deterministic Markov processes (PDMPs).
#'
#' @param nSample Number of samples after burn-in.
#' @param burnin Number of burn-in samples (default = 0).
#' @param mean A d-dimensional mean vector.
#' @param prec A d-by-d precision matrix of the Gaussian distribution.
#' @param lowerBounds A d-dimensional vector specifying the lower bounds.
#'   `-Inf` is accepted.
#' @param upperBounds A d-dimensional vector specifying the upper bounds.
#'   `Inf` is accepted.
#' @param init A d-dimensional vector of the initial value. `init` must
#'   satisfy all constraints. If `init = NULL`, a random initial value will
#'   be used.
#' @param stepSize Step size for the Markovian Zigzag sampler. Default value 
#'   is the empirically optimal choice: \eqn{\sqrt{2}\lambda^{-1/2}}, where \eqn{\lambda} is the 
#'   minimal eigenvalue of the precision matrix.
#' @param seed Random seed (default = 1).
#' @param numThreads number of threads for parallel execution (default = 1). Set to 0 for automatic detection of available cores.
#' @param diagnosticMode Logical. `TRUE` for also returning diagnostic
#'   information such as the step size used.
#' @param nStatusUpdate Number of status updates to print during sampling.
#'   If 0 (default), no updates are printed.
#'
#' @return An nSample-by-d matrix of samples. If `diagnosticMode` is `TRUE`,
#'   a list with additional diagnostic information is returned.
#' @export
#' @examples
#' set.seed(1)
#' d <- 5
#' A <- matrix(runif(d^2)*2-1, ncol=d)
#' precMat <- t(A) %*% A
#' initial <- rep(1, d)
#' results <- markovianZigzag(
#'   nSample = 1000,
#'   burnin = 1000,
#'   mean = rep(0, d),
#'   prec = precMat,
#'   lowerBounds = rep(0, d),
#'   upperBounds = rep(Inf, d)
#' )
#' @references
#' Bierkens, J., Roberts, G. O., and Zitt, P.-A. (2019). Ergodicity of the 
#' zigzag process. The Annals of Applied Probability, 29(4): 2266-2301.
#' @seealso [getMarkovianZigzagSample()], [createEngine()]

markovianZigzag <- function(nSample,
                            burnin = 0,
                            mean,
                            prec,
                            lowerBounds,
                            upperBounds,
                            init = NULL,
                            stepSize = NULL,
                            seed = 1,
                            numThreads = 1,
                            diagnosticMode = FALSE,
                            nStatusUpdate = 0L) {
  validateInput(mean, prec, lowerBounds, upperBounds, init)
  if (is.null(init)) {
    init <- getInitialPosition(mean, lowerBounds, upperBounds)
  }
  nIterPerUpdate <- ceiling((nSample + burnin) / nStatusUpdate)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  cpp_seed <- sample.int(.Machine$integer.max, size = 1)
  ndim <- length(mean)
  samples <- array(0, c(nSample, ndim))
  
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
  
  velocity <- 2 * stats::rbinom(ndim, 1, .5) - 1
  state <- list(position = init, velocity = velocity)
  for (i in 1:(nSample + burnin)) {
    state <- getMarkovianZigzagSample(
      position = state$position,
      velocity = state$velocity,
      engine = engine,
      travelTime = stepSize
    )
    if (i > burnin) {
      samples[i - burnin, ] <- state$position
    }
    if (i %% nIterPerUpdate == 0) {
      print(sprintf("%s iterations completed.", as.integer(i)))
    }
  }
  if (diagnosticMode) {
    return(list("samples" = samples, "stepsize" = stepSize))
  } else {
    return(samples)
  }
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
#' @examples
#' # First create an engine
#' set.seed(123)
#' engine <- createEngine(
#'   dimension = 2,
#'   lowerBounds = c(-1, -1),
#'   upperBounds = c(1, 1),
#'   seed = 123,
#'   mean = c(0, 0),
#'   precision = diag(2)
#' )
#' 
#' # Draw a single Markovian zigzag sample
#' position <- c(0.1, -0.2)
#' travel_time <- 0.5
#' sample_result <- getMarkovianZigzagSample(
#'   position = position,
#'   engine = engine,
#'   travelTime = travel_time
#' )
#' sample_result

#' @seealso [markovianZigzag()]
getMarkovianZigzagSample <- function(position, velocity = NULL, engine, travelTime) {
  if (is.null(velocity)) {
    velocity <- 2 * stats::rbinom(length(position), 1, .5) - 1
  }
  
  res <- .oneIrreversibleIteration(
    sexp = engine$engine,
    position = position,
    velocity = velocity,
    time = travelTime
  )
  return(list(
    position = res$position,
    velocity = res$velocity
  ))
}