#' Run Harmonic HMC Sampler
#'
#' Sample from a truncated Gaussian distribution with constraints Fx+g >= 0.
#'
#' @param n number of samples
#' @param initialPosition starting value for parameters
#' @param constraintDirec F matrix (k-by-d matrix where k is the number of
#' linear constraints)
#' @param constraintBound g vector (k dimensional)
#' @param choleskyFactor upper triangular matrix R from cholesky decomposition
#' of precision or covariance matrix into R^TR
#' @param unconstrainedMean mean of unconstrained Gaussian
#' @param precParametrized boolean for whether parametrization is by precision
#' (TRUE) or covariance matrix (FALSE)
#' @param integrationTime amount of time the particle travels for each sample.
#' Can either be a scalar value for a fixed time across all samples, or a length 2
#' vector of a lower and upper bound for uniform distribution from which the
#' bounce time is drawn from for each sample.
#' @param seed random seed
#' @param diagnosticMode boolean for whether to return the bounce distances for
#' each sample
#' @return List of
#' "samples": d x n matrix of samples
#' "bounceDistances": list of bounces for each sample, only present if
#' diagnosticMode is TRUE
#' @export
#'
#' @examples
#' set.seed(1)
#' d = 10
#' A = matrix(runif(d^2)*2-1, ncol=d)
#' Sigma = t(A) %*% A
#' R = cholesky(Sigma)
#' mu = rep(0,d)
#' constraintDirec = diag(d)
#' constraintBound = rep(0,d)
#' initial = rep(1, d)
#' results = runHHMC(
#' 100,
#' initial,
#' constraintDirec,
#' constraintBound,
#' R,
#' mu,
#' precParametrized = FALSE
#' )

runHHMC = function(n,
                   initialPosition,
                   constraintDirec,
                   constraintBound,
                   choleskyFactor,
                   unconstrainedMean,
                   precParametrized = TRUE,
                   integrationTime = c(pi / 8, pi / 2),
                   seed = 1,
                   diagnosticMode = FALSE) {
  if (length(integrationTime)==1){
    integrationTime[2] = integrationTime[1]
  }
  if (integrationTime[2] < integrationTime[1]){
    stop("Upper bound for integration time must be greater than lower bound.")
  }
  set.seed(seed)
  samples = matrix(nrow = ncol(constraintDirec), ncol = n)
  randomBounceTime = ifelse(length(integrationTime)==2, TRUE, FALSE)
  bounceDistances = vector(mode = "list",
                           length = ifelse(diagnosticMode, n, 0))
  whitenedConstraints = applyWhitenTransform(
    constraintDirec,
    constraintBound,
    choleskyFactor,
    unconstrainedMean,
    precParametrized
  )
  position = whitenPosition(
    initialPosition,
    constraintDirec,
    constraintBound,
    choleskyFactor,
    unconstrainedMean,
    precParametrized
  )
  for (i in 1:n) {
    momentum = rnorm(ncol(constraintDirec))
    results =  simulateWhitenedDynamics(
      position,
      momentum,
      whitenedConstraints$direc,
      whitenedConstraints$direcRowNormSq,
      whitenedConstraints$bound,
      runif(1, integrationTime[1], integrationTime[2]),
      diagnosticMode
    )
    position = results$position
    samples[, i] = unwhitenPosition(position,
                                    choleskyFactor,
                                    unconstrainedMean,
                                    precParametrized)
    if (diagnosticMode) {
      bounceDistances[[i]] = results$bounceDistances
    }
  }
  if (diagnosticMode) {
    return(list("samples" = samples, "bounceDistances" = bounceDistances))
  } else {
    return(list("samples" = samples))
  }
}
