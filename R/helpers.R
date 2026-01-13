#' Draw a random Laplace momentum
#'
#' Generate a d-dimensional momentum where the density of each element is proportional to exp(-|pi|).  
#' 
#' @param d dimension of the momentum.
#'
#' @return a d-dimensional Laplace-distributed momentum.
#' @export
#' @examples
#' # Draw a 3-dimensional Laplace momentum with reproducible results
#' set.seed(3)
#' momentum <- drawLaplaceMomentum(3)
#' momentum
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
#' @param upperBounds a d-dimensional vector specifying the upper bounds.
#'
#' @return an eligible d-dimensional initial vector.
#' @export
#' @examples
#' # Example 1: Bounded interval
#' mean <- c(0, 0)
#' lower <- c(-1, -2)
#' upper <- c(1, 2)
#' getInitialPosition(mean, lower, upper)
#'
#' # Example 2: Mixed bounds (some finite, some infinite)
#' mean <- c(0, 0, 0)
#' lower <- c(-Inf, 0, -1)
#' upper <- c(Inf, 5, Inf)
#' getInitialPosition(mean, lower, upper)
#'
#' # Example 3: All unbounded (returns mean)
#' mean <- c(1, 2, 3)
#' lower <- c(-Inf, -Inf, -Inf)
#' upper <- c(Inf, Inf, Inf)
#' getInitialPosition(mean, lower, upper)
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

#' Compute extreme eigenvalues of a symmetric matrix
#' 
#' Internal function to compute either the smallest or largest eigenvalue
#' of a symmetric matrix using Lanczos algorithm.
#'
#' @param symMatrix Symmetric matrix
#' @param smallest Logical, if TRUE computes smallest eigenvalue, 
#'        if FALSE computes largest eigenvalue
#' @param tol Tolerance for Lanczos algorithm convergence
#' @return The requested extreme eigenvalue
#' @noRd
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


#' Validate input parameters for MTN sampling
#' 
#' Internal validation function that checks consistency of parameters for 
#' multivariate truncated normal sampling.
#'
#' @param mean Mean vector
#' @param prec Precision matrix
#' @param lowerBounds Vector of lower bounds
#' @param upperBounds Vector of upper bounds  
#' @param init Optional initial position vector
#' @return Invisible NULL if validation passes, otherwise stops with error
#' @noRd
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
