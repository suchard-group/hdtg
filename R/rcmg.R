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
#' @param cppFlg
#'
#' @return An n-by-d matrix where each row is a multivariate sample
#' @export

rcmg <- function(n,
                 mean,
                 cov = NULL,
                 prec = NULL,
                 constraits,
                 t,
                 burnin,
                 p0 = NULL,
                 fixedMomentum = NULL,
                 cppFlg = FALSE,
                 nutsFlg = FALSE,
                 random_seed = 666,
                 randomFlg = TRUE,
                 debug_flg = F) {
  require(matrixcalc)
  stopifnot("n > burnin must be integers!" = (n %% 1 == 0 &&
                                                burnin %% 1 == 0 &&
                                                n > burnin))
  
  stopifnot("must provide mean" = !is.null(mean))
  
  
  if(!is.null(prec)){
    stopifnot("precision matrix contains NaN or is not positive definite " = (!any(is.na(prec)) && is.positive.definite(prec)))
  } else if (!is.null(prec)){
    stopifnot("covariance matrix contains NaN or is not positive definite " = (!any(is.na(cov)) && is.positive.definite(cov)))
    prec = solve(cov)
  } else {
    stop("must provide precision or covariance matrix")
  }
  
  if (!is.null(p0)) {
    stopifnot("initial value p0 does not comply with constraits" = all(p0 * constraits >= 0))
  } else {
    p0 = getInitialValue(mean, constraits)
  }
  # TODO add other checks for arguments. all dimensions must match.
  
  ndim = length(mean)
  energyGrad = function (x) {
    if (length(x) == 1) {
      return(prec[, x])
    } else {
      return(drop(prec %*% x))
    }
  }
  
  samples = array(0, c(n, ndim))
  set.seed(random_seed)
  
  if (cppFlg) {
    if (nutsFlg) {
      engine = createNutsEngine(
        dimension = ndim,
        mask = rep(1, ndim),
        observed = rep(1, ndim),
        parameterSign = constraits,
        flags = 128,
        info = 1,
        seed = random_seed,
        randomFlg = randomFlg,
        stepSize = t,
        mean = mean,
        precision = prec
      )
    } else {
      engine = createEngine(
        dimension = ndim,
        mask = rep(1, ndim),
        observed = rep(1, ndim),
        parameterSign = constraits,
        flags = 128,
        info = 1,
        seed = random_seed)
    }
  }
  
  if (cppFlg) {
    for (i in 1:n) {
      if (!is.null(fixedMomentum)) {
        momentum = fixedMomentum
      } else {
        momentum = drawMomentum(ndim)
      }
      
      p0 = getSample(
        position = p0,
        momentum = momentum,
        t = t,
        nutsFlg = nutsFlg,
        engine = engine
      )
      samples[i, ] = p0
      
      if (debug_flg) {
        cat('iteration', i, 'done \n')
      }
    }
  } else {
    for (i in 1:n) {
      if (!is.null(fixedMomentum)) {
        momentum = fixedMomentum
      } else {
        momentum = drawMomentum(ndim)
      }
      
      p0 = getSampleR(
        position = p0,
        momentum = momentum,
        t = t,
        energyGrad = energyGrad,
        mean = mean,
        constraits = constraits
      )
      
      samples[i, ] = p0
      
      if (debug_flg) {
        cat('iteration', i, 'done \n')
      }
    }
  }
  return(samples)
}

#' @param mean
#'
#' @param constraits
#'
#' @export
getInitialValue <- function(mean, constraits) {
  p0 = mean
  p0[p0 == 0] = 0.1 * constraits[p0 == 0]
  p0[p0 * constraits < 0] = constraits[p0 * constraits < 0]
  return(p0)
}

#' Title
#'
#' @param dim 
#'
#' @return
#' @export
#'
#' @examples
drawMomentum <- function(dim){
  return((2 * (runif(dim) > .5) - 1) * rexp(dim, rate = 1))
}
