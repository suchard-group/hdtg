#' Sample from a truncated Gaussian distribution with the harmonic HMC
#'
#' Generate MCMC samples from a d-dimensional truncated Gaussian distribution
#' with constraints Fx+g >= 0 using the Harmonic Hamiltonian Monte Carlo sampler 
#' (Harmonic-HMC).
#' @param n number of samples after burn-in.
#' @param burnin number of burn-in samples.
#' @param mean a d-dimensional mean vector.
#' @param choleskyFactor upper triangular matrix R from Cholesky decomposition
#' of precision or covariance matrix into R^TR.
#' @param constraintDirec F matrix (k-by-d matrix where k is the number of
#' linear constraints).
#' @param constraintBound g vector (k-dimensional).
#' @param init a d-dimensional vector of the initial value. `init` must satisfy
#' all constraints.
#' @param integrationTime HMC integration time for each iteration. Can either be
#' a scalar value for a fixed time across all samples, or a length 2 vector of a
#' lower and upper bound for uniform distribution from which the time is drawn
#' from for each iteration.
#' @param precParametrized logical. whether `choleskyFactor` is from precision
#' (`TRUE`) or covariance matrix (`FALSE`).
#' @param diagnosticMode logical. `TRUE` for also returning the bounce distances
#' for each sample.
#' @return List of
#' `samples`: (n + burnin) x d matrix of samples (including burnin samples) and
#' `bounceDistances`: list of bounces for each sample (only present if
#' `diagnosticMode` is `TRUE`).
#' @export
#'
#' @examples
#' set.seed(1)
#' d <- 10
#' A <- matrix(runif(d^2)*2-1, ncol=d)
#' Sigma <- t(A) %*% A
#' R <- cholesky(Sigma)
#' mu <- rep(0,d)
#' constraintDirec <- diag(d)
#' constraintBound <- rep(0,d)
#' initial <- rep(1, d)
#' results <- harmonicHMC(
#' 1000,
#' 1000,
#' mu,
#' R,
#' constraintDirec,
#' constraintBound,
#' initial,
#' precParametrized = FALSE)
#' @references
#' \insertRef{pakman2014exact}{hdtg}

harmonicHMC <- function(n,
                        burnin,
                        mean,
                        choleskyFactor,
                        constraintDirec,
                        constraintBound,
                        init,
                        integrationTime = c(pi / 8, pi / 2),
                        precParametrized = TRUE,
                        diagnosticMode = FALSE) {
  if (length(integrationTime) == 1) {
    integrationTime[2] <- integrationTime[1]
  }
  if (integrationTime[2] < integrationTime[1]) {
    stop("Upper bound for integration time must be greater than lower bound.")
  }
  if (sum(constraintDirec %*% init + constraintBound) < length(constraintBound)) {
    stop("Initial value x does not satisfy Fx + g >=0")
  }
  samples <- matrix(nrow = n + burnin, ncol = ncol(constraintDirec))
  randomBounceTime <-
    ifelse(length(integrationTime) == 2, TRUE, FALSE)
  bounceDistances <- vector(mode = "list",
                            length = ifelse(diagnosticMode, n + burnin, 0))
  whitenedConstraints <- applyWhitenTransform(constraintDirec,
                                              constraintBound,
                                              choleskyFactor,
                                              mean,
                                              precParametrized)
  position <- whitenPosition(init,
                             constraintDirec,
                             constraintBound,
                             choleskyFactor,
                             mean,
                             precParametrized)
  for (i in 1:(n + burnin)) {
    momentum <- rnorm(ncol(constraintDirec))
    results <- simulateWhitenedDynamics(
      position,
      momentum,
      whitenedConstraints$direc,
      whitenedConstraints$direcRowNormSq,
      whitenedConstraints$bound,
      runif(1, integrationTime[1], integrationTime[2]),
      diagnosticMode
    )
    position <- results$position
    samples[i,] <- unwhitenPosition(position,
                                    choleskyFactor,
                                    mean,
                                    precParametrized)
    if (diagnosticMode) {
      bounceDistances[[i]] <- results$bounceDistances
    }
  }
  if (diagnosticMode) {
    return(list("samples" = samples, "bounceDistances" = bounceDistances))
  } else {
    return(list("samples" = samples))
  }
}
