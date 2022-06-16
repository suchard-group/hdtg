#' Sample from a truncated Gaussian distribution with the zigzag HMC
#'
#' Generate MCMC samples from a d-dimensional truncated Gaussian distribution with element-wise truncations using the Zigzag Hamiltonian Monte Carlo sampler (Zigzag-HMC).
#'
#' @param n number of samples after burn-in.
#' @param burnin number of burn-in samples.
#' @param mean a d-dimensional mean vector.
#' @param cov a d-by-d covariance matrix of the Gaussian distribution. At least one of `prec` and `cov` should be provided.
#' @param prec a d-by-d precision matrix of the Gaussian distribution. 
#' @param lowerBounds a d-dimensional vector specifying the lower bounds. `-Inf` is accepted.  
#' @param upperBounds a d-dimensional vector specifying the upper bounds. `Inf` is accepted. 
#' @param nutsFlg logical. If `TRUE` the No-U-Turn sampler will be used (Zigzag-NUTS).
#' @param init a d-dimensional vector of the initial value. `init` must satisfy all constraints. If `init = NULL`, a random initial value will be used.
#' @param step step size for Zigzag-HMC or Zigzag-NUTS (if `nutsFlg = TRUE`). Default value is the emipirically optimal choice: sqrt(2)(lambda)^(-1/2) for Zigzag-HMC and 0.1(lambda)^(-1/2) for Zigzag-NUTS, where lambda is the minimal eigenvalue of the precision matrix.   
#' @param rSeed random seed (default = 1).
#'
#' @return (n + burnin) x d matrix of samples. The first `burnin` samples are from the user specified warm-up iterations.
#' @export
#' 
#' @examples
#' set.seed(1)
#' d <- 10
#' A <- matrix(runif(d^2)*2-1, ncol=d)
#' covMat <- t(A) %*% A
#' initial <- rep(1, d)
#' results <- zigzagHMC(
#' n = 1000,
#' burnin = 1000,
#' mean = rep(0, d),
#' cov = covMat,
#' lowerBounds = rep(0, d),
#' upperBounds = rep(Inf, d))
#'
#' @references
#' \insertRef{nishimura2021hamiltonian}{hdtg}
#'
#' \insertRef{nishimura2020discontinuous}{hdtg}

zigzagHMC <- function(n,
                      burnin,
                      mean,
                      cov,
                      prec = NULL,
                      lowerBounds,
                      upperBounds,
                      init = NULL,
                      step = NULL,
                      nutsFlg = FALSE,
                      rSeed = 1) {
  ndim <- length(mean)
  
  if (!is.null(prec)) {
    stopifnot("precision matrix contains NaN" = !any(is.na(prec)))
  } else if (!is.null(cov)) {
    stopifnot("covariance matrix contains NaN" = !any(is.na(cov)))
    prec <- solve(cov)
  } else {
    stop("must provide precision or covariance matrix")
  }
  
  stopifnot(
    "precision / covariance matrix has incompatible dimensions" = (nrow(prec) == ndim &&
                                                                     ncol(prec) == ndim)
  )
  stopifnot(
    "some lower bound is larger than the corresponding upper bound" = sum(lowerBounds < upperBounds) == ndim
  )
  
  if (!is.null(init)) {
    stopifnot(
      "initial position is not compatiable with the truncation bounds" = (sum(lowerBounds < init) == ndim) &&
        (sum(init < upperBounds) == ndim)
    )
  } else {
    init <- getInitialPosition(mean, lowerBounds, upperBounds)
  }
  
  energyGrad <- function (x) {
    if (length(x) == 1) {
      return(prec[, x])
    } else {
      return(drop(prec %*% x))
    }
  }
  
  samples <- array(0, c(n + burnin, ndim))
  
  if (nutsFlg) {
    if (!is.null(step)) {
      t <- step
    } else{
      t <- 0.1 / sqrt(min(mgcv::slanczos(
        A = prec, k = 1, kl = 1
      )[['values']]))
    }
    cat("NUTS base step size is", t, "\n")
    engine <- createNutsEngine(
      dimension = ndim,
      lowerBounds = lowerBounds,
      upperBounds = upperBounds,
      flags = 128,
      info = 1,
      seed = rSeed,
      randomFlg = FALSE,
      stepSize = t,
      mean = mean,
      precision = prec
    )
    
  } else {
    if (!is.null(step)) {
      t <- step
    } else{
      t <-
        sqrt(2) / sqrt(min(mgcv::slanczos(
          A = prec, k = 1, kl = 1
        )[['values']], na.rm = T))
    }
    cat("HZZ step size is", t, "\n")
    engine <- createEngine(
      dimension = ndim,
      lowerBounds = lowerBounds,
      upperBounds = upperBounds,
      flags = 128,
      info = 1,
      seed = rSeed
    )
    setMean(sexp = engine$engine, mean = mean)
    setPrecision(sexp = engine$engine, precision = prec)
  }
  
  position <- init
  for (i in 1:(n + burnin)) {
    position <- getSample(
      position = position,
      momentum = NULL,
      t = t,
      nutsFlg = nutsFlg,
      engine = engine
    )
    samples[i, ] <- position
  }
  return(samples)
}



#' A function to get an eligible initial value for a MTN
#'
#' @param mean a d-dimensional mean vector
#' @param lowerBounds a d-dimensional lower bound
#' @param upperBounds a d-dimensional lower bound
#'
#' @return an eligible d-dimensional initial value
#' @export
#'
getInitialPosition <- function(mean, lowerBounds, upperBounds) {
  bL <- upperBounds - lowerBounds
  midPoint <- (upperBounds + lowerBounds) / 2
  x <- mean
  x[is.finite(bL)] = midPoint[is.finite(bL)]
  x[is.infinite(bL) &
      is.finite(lowerBounds)] = lowerBounds[is.infinite(bL) &
                                              is.finite(lowerBounds)] + 0.1
  x[is.infinite(bL) &
      is.finite(upperBounds)] = upperBounds[is.infinite(bL) &
                                              is.finite(upperBounds)] - 0.1
  return(x)
}
