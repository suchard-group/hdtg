#' Sample from a truncated Gaussian distribution with the harmonic HMC
#'
#' Generate MCMC samples from a d-dimensional truncated Gaussian distribution
#' with constraints Fx+g >= 0 using the Harmonic Hamiltonian Monte Carlo sampler
#' (Harmonic-HMC).
#'
#' @param nSample number of samples after burn-in.
#' @param burnin number of burn-in samples (default = 0).
#' @param mean a d-dimensional mean vector.
#' @param choleskyFactor upper triangular matrix R from Cholesky decomposition of
#' precision or covariance matrix into R^TR.
#' @param constrainDirec the k-by-d F matrix (k is the number of linear constraints).
#' @param constrainBound the k-dimensional g vector.
#' @param init a d-dimensional vector of the initial value. `init` must satisfy all constraints.
#' @param time HMC integration time for each iteration. Can either be
#' a scalar value for a fixed time across all samples, or a length 2 vector of a
#' lower and upper bound for uniform distribution from which the time is drawn
#' from for each iteration.
#' @param precFlg logical. whether `choleskyFactor` is from precision
#' (`TRUE`) or covariance matrix (`FALSE`).
#' @param seed random seed (default = 1).
#' @param extraOutputs vector of strings. "numBounces" and/or "bounceDistances"
#' can be requested, with the latter containing the distances in-between bounces
#' for each sample and hence incurring significant computational and memory costs.
#' @return 
#' When `extraOutputs` is empty (default), returns an `nSample`-by-`d` matrix of samples.
#' 
#' When `extraOutputs` contains `"numBounces"` and/or `"bounceDistances"`, returns a list with elements:
#' \item{samples}{`nSample`-by-`d` matrix of samples}
#' \item{numBounces}{Vector of bounce counts per sample (if requested)}
#' \item{bounceDistances}{List of bounce distances per sample (if requested)}
#' @export
#'
#' @examples
#' set.seed(1)
#' d <- 10
#' A <- matrix(runif(d^2)*2 - 1, ncol=d)
#' precMat <- t(A) %*% A
#' R <- cholesky(precMat)
#' mu <- rep(0, d)
#' constrainDirec <- diag(d)
#' constrainBound <- rep(0,d)
#' initial <- rep(1, d)
#' results <- harmonicHMC(1000, 1000, mu, R, constrainDirec, constrainBound, initial, precFlg = TRUE)
#' @references
#'
#' Pakman, A. and Paninski, L. (2014). Exact Hamiltonian Monte Carlo for 
#' Truncated Multivariate Gaussians. Journal of Computational and Graphical 
#' Statistics. doi:10.1080/10618600.2013.788448
#' 
#' @seealso [getHarmonicSample()], [cholesky()], [getInitialPosition()]
harmonicHMC <- function(nSample,
                        burnin = 0,
                        mean,
                        choleskyFactor,
                        constrainDirec,
                        constrainBound,
                        init,
                        time = c(pi / 8, pi / 2),
                        precFlg,
                        seed = 1,
                        extraOutputs = c()) {
  if (length(time) == 1) {
    time[2] <- time[1]
  }
  if (time[2] < time[1]) {
    stop("Upper bound for integration time must be greater than lower bound.")
  }
  if (sum(constrainDirec %*% init + constrainBound > 0) < length(constrainBound)) {
    stop("Initial value x does not satisfy Fx + g >=0")
  }
  samples <- matrix(nrow = nSample, ncol = ncol(constrainDirec))
  randomBounceTime <-
    ifelse(length(time) == 2, TRUE, FALSE)
  numBounces <- vector(mode = "integer", length = nSample)
  diagnosticMode <- "bounceDistances" %in% extraOutputs
  bounceDistances <- vector(mode = "list",
                            length = ifelse(diagnosticMode, nSample, 0))
  whitenedConstraints <- applyWhitenTransform(constrainDirec,
                                              constrainBound,
                                              choleskyFactor,
                                              mean,
                                              precFlg)
  position <- whitenPosition(init,
                             constrainDirec,
                             constrainBound,
                             choleskyFactor,
                             mean,
                             precFlg)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  for (i in 1:(nSample + burnin)) {

    results <- getHarmonicSample(
      whitenedPosition = position,
      whitenedConstraints = whitenedConstraints,
      integrationTime = runif(1, time[1], time[2]),
      diagnosticMode = diagnosticMode
    )
    
    position <- results$position
    
    if (i > burnin) {
      samples[i - burnin, ] <- unwhitenPosition(position, choleskyFactor, mean, precFlg)
      numBounces[i - burnin] <- results$numBounces
      if ("bounceDistances" %in% extraOutputs) {
        bounceDistances[[i - burnin]] <- results$bounceDistances
      }
    }
  }
  if (length(extraOutputs) == 0) {
    return(samples)
  } else {
    output <- list("samples" = samples)
    for (quantity in extraOutputs) {
      output[[quantity]] <- get(quantity)
    }
    return(output)
  }
}


#' One-step Harmonic HMC Sampler (Whitened Coordinates)
#'
#' @param whitenedPosition Position in whitened coordinates
#' @param whitenedConstraints List from \code{applyWhitenTransform()}
#' @param integrationTime Time for dynamics simulation
#' @param diagnosticMode Return bounce diagnostics
#' @export
#' @examples
#' # Basic usage with whitened coordinates
#' set.seed(123)
#' whitened_pos <- c(0.1, -0.2, 0.3)
#' # Create example whitened constraints
#' whitened_constraints <- list(
#'   direc = matrix(c(1, 0, 0, 0, 1, 0), nrow = 2, byrow = TRUE),
#'   direcRowNormSq = c(1, 1),
#'   bound = c(-0.5, -0.5)
#' )
#' result <- getHarmonicSample(
#'   whitenedPosition = whitened_pos,
#'   whitenedConstraints = whitened_constraints,
#'   integrationTime = pi/4
#' )
#' result
#' 
#' # With diagnostics enabled
#' result_diag <- getHarmonicSample(
#'   whitenedPosition = whitened_pos,
#'   whitenedConstraints = whitened_constraints,
#'   integrationTime = pi/4,
#'   diagnosticMode = TRUE
#' )
#' str(result_diag)
#' @seealso [harmonicHMC()]
getHarmonicSample <- function(whitenedPosition,
                              whitenedConstraints,
                              integrationTime,
                              diagnosticMode = FALSE) {
  momentum <- rnorm(length(whitenedPosition))
  
  simulateWhitenedDynamics(
    whitenedPosition,
    momentum,
    whitenedConstraints$direc,
    whitenedConstraints$direcRowNormSq,
    whitenedConstraints$bound,
    integrationTime,
    diagnosticMode
  )
}