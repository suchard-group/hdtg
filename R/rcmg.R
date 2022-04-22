#' Sample from a constrained Gaussian distribution
#'
#' Generate samples from a constrained Gaussian distribution with the Hamiltonian zigzag sampler.
#'
#' @param n number of samples
#' @param mean a d-dimensional mean vector
#' @param prec a d-by-d precision matrix of the Gaussian distribution
#' @param constraits a d-dimensional vector where a positive (negative) number means >0 (<0) truncation.
#' @param p0 a d-dimensional vector of the initial value. It must satisfy all constraints. If not specified a random initial value will be used
#' @param t time length to simulate the Markov process
#' @param cppFlg
#'
#' @return An n-by-d matrix where each row is a multivariate sample
#' @export

rcmg <- function(n,
                 mean,
                 cov,
                 prec = NULL,
                 lowerBounds,
                 upperBounds,
                 #burnin,
                 p0 = NULL,
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
  
  # if (!is.null(p0)) {
  #   stopifnot("initial value p0 does not comply with constraits" = all(p0 * constraits >= 0))
  # } else {
  #   p0 <- getInitialValueNew(mean, lowerBounds, upperBounds)
  # }
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
    if(!is.null(forcedStep)){
      t <- forcedStep
    }else{
      t <- 0.1 / sqrt(min(mgcv::slanczos(A = prec,k=1,kl=1)[['values']]))
    }
    cat("NUTS base step size is", t)
    engine <- createNutsEngine(
      dimension = ndim,
      mask = rep(1, ndim),
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
    if(!is.null(forcedStep)){
      t <- forcedStep
    }else{
      t <- sqrt(2) / sqrt(min(mgcv::slanczos(A = prec,k=1,kl=1)[['values']], na.rm = T))
    }
    cat("HZZ step size is", t)
    engine <- createEngine(
      dimension = ndim,
      mask = rep(1, ndim),
      lowerBounds = lowerBounds,
      upperBounds = upperBounds,
      flags = 128,
      info = 1,
      seed = rSeed)
    setMean(sexp = engine$engine, mean = mean)
    setPrecision(sexp = engine$engine, precision = prec)
  }
  
  
  set.seed(rSeed)
  
  for (i in 1:n) {
    p0 <- getSample(
      position = p0,
      momentum = momentum,
      t = t,
      nutsFlg = nutsFlg,
      engine = engine
    )
    samples[i, ] <- p0
  }
  return(samples)
}

#' @param mean
#'
#' @param constraits
#'
#' @export
getInitialValue <- function(mean, constraits) {
  p0 <- mean
  p0[p0 == 0] <- 0.1 * constraits[p0 == 0]
  p0[p0 * constraits < 0] <- constraits[p0 * constraits < 0]
  return(p0)
}

#' @param mean the mean vector for the MTN distribution
#'
#' @param lowerBounds the lower bound vector
#' @param upperBounds the upper bound vector
#'
#' @export
getInitialValueNew <- function(mean, lowerBounds, upperBounds) {
  p0 <- mean
  p0[mean >= upperBounds] <- upperBounds[mean >= upperBounds]
  p0[mean <= lowerBounds] <- lowerBounds[mean <= upperBounds]
  return(p0)
}

