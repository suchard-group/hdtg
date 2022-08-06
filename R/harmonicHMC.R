#' Sample from a truncated Gaussian distribution with the harmonic HMC
#'
#' Generate MCMC samples from a d-dimensional truncated Gaussian distribution
#' with constraints Fx+g >= 0 using the Harmonic Hamiltonian Monte Carlo sampler 
#' (Harmonic-HMC).
#'
#' @param n number of samples after burn-in.
#' @param burnin number of burn-in samples (default = 0).
#' @param mean a d-dimensional mean vector.
#' @param choleskyFactor upper triangular matrix R from Cholesky decomposition of 
#' precision or covariance matrix into R^TR. 
#' @param F F matrix (k-by-d matrix where k is the number of linear constraints). 
#' @param g g vector (k-dimensional). 
#' @param init a d-dimensional vector of the initial value. `init` must satisfy all constraints. 
#' @param time HMC integration time for each iteration. Can either be
#' a scalar value for a fixed time across all samples, or a length 2 vector of a
#' lower and upper bound for uniform distribution from which the time is drawn
#' from for each iteration.
#' @param precFlg logical. whether `choleskyFactor` is from precision
#' (`TRUE`) or covariance matrix (`FALSE`).
#' @param diagnosticMode logical. `TRUE` for also returning the bounce distances
#' for each sample.
#'
#' @return List of
#' `samples`: (n + burnin) x d matrix of samples (including burnin samples) and
#' `bounceDistances`: list of bounces for each sample (only present if
#' `diagnosticMode` is `TRUE`).
#' @export
#'
#' @examples
#' set.seed(1)
#' d <- 10
#' A <- matrix(runif(d^2)*2 - 1, ncol=d)
#' Sigma <- t(A) %*% A
#' R <- cholesky(Sigma)
#' mu <- rep(0, d)
#' F <- diag(d)
#' g <- rep(0,d)
#' initial <- rep(1, d)
#' results <- harmonicHMC(1000, 1000, mu, R, F, g, initial, precFlg = FALSE)
#' @references
#' \insertRef{pakman2014exact}{hdtg}
harmonicHMC <- function(n,
                        burnin = 0,
                        mean,
                        choleskyFactor,
                        F,
                        g,
                        init,
                        time = c(pi / 8, pi / 2),
                        precFlg,
                        diagnosticMode = FALSE) {
  if (length(time) == 1) {
    time[2] <- time[1]
  }
  if (time[2] < time[1]) {
    stop("Upper bound for integration time must be greater than lower bound.")
  }
  if (sum(F %*% init + g > 0) < length(g)) {
    stop("Initial value x does not satisfy Fx + g >=0")
  }
  samples <- matrix(nrow = n + burnin, ncol = ncol(F))
  randomBounceTime <-
    ifelse(length(time) == 2, TRUE, FALSE)
  bounceDistances <- vector(mode = "list",
                            length = ifelse(diagnosticMode, n + burnin, 0))
  whitenedConstraints <- applyWhitenTransform(F,
                                              g,
                                              choleskyFactor,
                                              mean,
                                              precFlg)
  position <- whitenPosition(init,
                             F,
                             g,
                             choleskyFactor,
                             mean,
                             precFlg)
  for (i in 1:(n + burnin)) {
    momentum <- rnorm(ncol(F))
    results <- simulateWhitenedDynamics(
      position,
      momentum,
      whitenedConstraints$direc,
      whitenedConstraints$direcRowNormSq,
      whitenedConstraints$bound,
      runif(1, time[1], time[2]),
      diagnosticMode
    )
    position <- results$position
    samples[i,] <- unwhitenPosition(position,
                                    choleskyFactor,
                                    mean,
                                    precFlg)
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
