#' Sample from a constrained Gaussian distribution
#'
#' Generate samples from a constrained Gaussian distribution with the Hamiltonian zigzag sampler.
#'
#' @param n number of samples
#' @param mean a d-dimensional mean vector
#' @param prec a d-by-d precision matrix of the Gaussian distribution
#' @param cov todo
#' @param lowerBounds todo
#' @param upperBounds todo
#' @param x a d-dimensional vector of the initial value. It must satisfy all constraints. If not specified a random initial value will be used
#' @param forcedStep todo
#' @param momentum todo
#' @param nutsFlg todo
#' @param rSeed todo
#' @param randomFlg todo
#'
#' @return An n-by-d matrix where each row is a multivariate sample
#' @export

rcmg <- function(n,
                 mean,
                 cov,
                 prec = NULL,
                 lowerBounds,
                 upperBounds,
                 x = NULL,
                 forcedStep = NULL,
                 momentum = NULL,
                 nutsFlg = FALSE,
                 rSeed = 666,
                 randomFlg = TRUE) {
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
  if (!is.null(x)) {
    stopifnot(
      "initial position is not compatiable with the truncation bounds" = (sum(lowerBounds < x) == ndim) &&
        (sum(x < upperBounds) == ndim)
    )
  } else {
    x <- getInitialPosition(mean, lowerBounds, upperBounds)
  }
  
  # TODO add other checks for arguments. all dimensions must match.
  
  energyGrad <- function (x) {
    if (length(x) == 1) {
      return(prec[, x])
    } else {
      return(drop(prec %*% x))
    }
  }
  
  samples <- array(0, c(n, ndim))
  
  if (nutsFlg) {
    if (!is.null(forcedStep)) {
      t <- forcedStep
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
      randomFlg = randomFlg,
      stepSize = t,
      mean = mean,
      precision = prec
    )
    
  } else {
    if (!is.null(forcedStep)) {
      t <- forcedStep
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
  
  for (i in 1:n) {
    x <- getSample(
      position = x,
      momentum = momentum,
      t = t,
      nutsFlg = nutsFlg,
      engine = engine
    )
    samples[i,] <- x
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
  x[is.infinite(bL) & is.finite(lowerBounds)] = lowerBounds[is.infinite(bL) & is.finite(lowerBounds)] + 0.1
  x[is.infinite(bL) & is.finite(upperBounds)] = upperBounds[is.infinite(bL) & is.finite(upperBounds)] - 0.1
  return(x)
}
