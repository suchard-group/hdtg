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
#' @param init a d-dimensional vector of the initial value. `init` must satisfy all constraints. If `init = NULL`, a random initial value will be used.
#' @param stepsize step size for Zigzag-HMC or Zigzag-NUTS (if `nutsFlg = TRUE`). Default value is the empirically optimal choice: sqrt(2)(lambda)^(-1/2) for Zigzag-HMC and 0.1(lambda)^(-1/2) for Zigzag-NUTS, where lambda is the minimal eigenvalue of the precision matrix.   
#' @param seed random seed (default = 1).
#' @param diagnosticMode logical. `TRUE` for also returning diagnostic information such as the stepsize used. 
#'
#' @return an nSample-by-d matrix of samples. If `diagnosticMode` is `TRUE`, a list with additional diagnostic information is returned. 
#' @export
#' @examples
#' set.seed(1)
#' d <- 10
#' A <- matrix(runif(d^2)*2-1, ncol=d)
#' covMat <- t(A) %*% A
#' initial <- rep(1, d)
#' results <- zigzagHMC(nSample = 1000, burnin = 1000, mean = rep(0, d), cov = covMat,
#' lowerBounds = rep(0, d), upperBounds = rep(Inf, d))
#'
#' @references
#' \insertRef{nishimura2021hamiltonian}{hdtg}
#'
#' \insertRef{nishimura2020discontinuous}{hdtg}

zigzagHMC <- function(nSample,
                      burnin = 0,
                      mean,
                      prec,
                      lowerBounds,
                      upperBounds,
                      init = NULL,
                      stepsize = NULL,
                      nutsFlg = FALSE,
                      seed = 1,
                      diagnosticMode = FALSE) {
  
  validateInput(mean, prec, lowerBounds, upperBounds, init)
  if (is.null(init)) {
    init <- getInitialPosition(mean, lowerBounds, upperBounds)
  }
  
  set.seed(seed)
  ndim <- length(mean)
  samples <- array(0, c(nSample, ndim))
  
  if (nutsFlg) {
    if (is.null(stepsize)) {
      stepsize <- 0.1 / sqrt(computeExtremeEigenval(Prec))
    }
    engine <- createNutsEngine(
      dimension = ndim,
      lowerBounds = lowerBounds,
      upperBounds = upperBounds,
      flags = 128,
      seed = seed,
      stepSize = stepsize,
      mean = mean,
      precision = prec
    )
    
  } else {
    if (is.null(stepsize)) {
      stepsize <- sqrt(2) / sqrt(computeExtremeEigenval(Prec))
    }
    engine <- createEngine(
      dimension = ndim,
      lowerBounds = lowerBounds,
      upperBounds = upperBounds,
      flags = 128,
      seed = seed,
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
      stepZZHMC = stepsize
    )
    if (i > burnin) {
      samples[i - burnin, ] <- position
    }
  }
  if (diagnosticMode) {
    return(list("samples" = samples, "stepsize" = stepsize))
  } else {
    return(samples)
  }
}

markovianZigzag <- function(nSample,
                            burnin = 0,
                            mean,
                            prec,
                            lowerBounds,
                            upperBounds,
                            init = NULL,
                            stepsize = NULL,
                            seed = 1,
                            diagnosticMode = FALSE,
                            nStatusUpdate = 0L) {
  
  validateInput(mean, prec, lowerBounds, upperBounds, init)
  if (is.null(init)) {
    init <- getInitialPosition(mean, lowerBounds, upperBounds)
  }
  nIterPerUpdate <- ceiling((nSample + burnin) / nStatusUpdate)
  
  set.seed(seed)
  ndim <- length(mean)
  samples <- array(0, c(nSample, ndim))
  
  if (is.null(stepsize)) {
    stepsize <- sqrt(2) / sqrt(computeExtremeEigenval(Prec))
  }
  engine <- createEngine(
    dimension = ndim,
    lowerBounds = lowerBounds,
    upperBounds = upperBounds,
    flags = 128,
    seed = seed,
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
      travelTime = stepsize
    )
    if (i > burnin) {
      samples[i - burnin, ] <- state$position
    }
    if (i %% nIterPerUpdate == 0) {
      print(sprintf("%s iterations completed.", as.integer(i)))
    }
  }
  if (diagnosticMode) {
    return(list("samples" = samples, "stepsize" = stepsize))
  } else {
    return(samples)
  }
}

computeExtremeEigenval <- function(symMatrix, smallest = TRUE, tol = .Machine$double.eps^.5) {
  if (smallest) {
    nLargest <- 0
    nSmallest <- 1
  } else {
    nLargest <- 1
    nSmallest <- 0
  }
  return(
    mgcv::slanczos(A = symMatrix, k = nLargest, kl = nSmallest, tol = tol)[['values']]
  )
}

validateInput <- function(mean, prec, lowerBounds, upperBounds, init) {
  ndim <- length(mean)
  stopifnot(
    "precision/covariance matrix size does not match the mean vector" = 
      (nrow(prec) == ndim && ncol(prec) == ndim)
  )
  stopifnot(
    "some lower bound is larger than the corresponding upper bound" = sum(lowerBounds < upperBounds) == ndim
  )
  if (!is.null(init)) {
    stopifnot(
      "initial position is not compatiable with the truncation bounds" = (sum(lowerBounds < init) == ndim) &&
        (sum(init < upperBounds) == ndim)
    )
  }
}
