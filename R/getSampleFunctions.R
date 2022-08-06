#' Simulate one MTN sample
#'  
#' Generate one random sample from the target MTN with Zigzag-HMC or Zigzag-NUTS.
#'
#' @param position a d-dimensional initial position vector.
#' @param momentum a d-dimensional initial momentum vector.
#' @param nutsFlg logical. If `TRUE` the No-U-Turn sampler will be used (Zigzag-NUTS).
#' @param engine list. Its `engine` element is a pointer to the Zigzag-HMC engine
#' (or Zigzag-NUTS engine) C++ object that implements fast computations for
#' Zigzag-HMC (or Zigzag-NUTS).
#' @param stepZZHMC step size for Zigzag-HMC. If `nutsFlg = TRUE`, `engine` contains
#' the base step size for Zigzag-NUTS).
#'
#' @return one MCMC sample from the target MTN.
#' @export
getZigzagSample <- function(position,
                            momentum = NULL,
                            nutsFlg,
                            engine,
                            stepZZHMC = NULL) {
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
      time = stepZZHMC
    )
  }
  return(res$position)
}

#' Draw a random Laplace momentum
#'
#' Generate a d-dimensional momentum where the density of each element is proportional to exp(-|pi|).  
#'
#' @param d dimension of the momentum.
#'
#' @return a d-dimensional Laplace-distributed momentum.
drawLaplaceMomentum <- function(d) {
  return((2 * (stats::runif(d) > .5) - 1) * stats::rexp(d, rate = 1))
}


#' Get an eligible initial value for a MTN with given mean and truncations
#'
#' For a given MTN the function returns an initial vector whose elements are one of:
#' (1) middle point of the truncation interval if both lower and upper bounds are
#' finite (2) lower (upper) bound +0.1 (-0.1) if only the lower (upper) bound is finite
#' (3) the corresponding mean value if lower bound = `-Inf` are upper bound = `Inf`.
#'
#' @param mean a d-dimensional mean vector.
#' @param lowerBounds a d-dimensional vector specifying the lower bounds.
#' @param upperBounds a d-dimensional vector specifying the lower bounds.
#'
#' @return an eligible d-dimensional initial vector.
#' @export
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