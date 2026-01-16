test_that("harmonicHMC matches TruncatedNormal reference for truncated Gaussian", {
  skip_if_not_installed("TruncatedNormal")
  library(TruncatedNormal)
  set.seed(123)
  
  # Define problem
  d <- 2
  meanVec <- c(0, 0)
  covMat <- matrix(c(1, 0.5, 0.5, 1), nrow = d)
  lb <- c(-1, -1)
  ub <- c(Inf, Inf)
  
  nSamples <- 50000
  burnin <- 20000
  nRef <- 100000
  
  R <- chol(covMat)
  constrainDirec <- diag(d)  
  constrainBound <- -lb  
  
  # Initial point that satisfies constraints
  init <- c(0.5, 0.5)
  
  # Reference samples
  samples_ref <- rtmvnorm(
    n = nRef,
    mu = meanVec,
    sigma = covMat,
    lb = lb,
    ub = ub
  )
  
  # harmonicHMC samples
  samples_hmc <- harmonicHMC(
    nSample = nSamples,
    burnin = burnin,
    mean = meanVec,
    choleskyFactor = R,
    constrainDirec = constrainDirec,
    constrainBound = constrainBound,
    init = init,
    precFlg = FALSE,  
    seed = 123,
    extraOutputs = c()
  )
  
  # Compare means
  ref_means <- colMeans(samples_ref)
  hmc_means <- colMeans(samples_hmc)
  
  # Print for debugging
  cat("\nReference means:", ref_means)
  cat("\nHarmonic HMC means:", hmc_means, "\n")
  
  tol <- 0.05
  expect_equal(hmc_means, ref_means, tolerance = tol)
})

test_that("harmonicHMC produces consistent results with seeding", {
  d <- 2
  mean_vec <- c(1, 2)
  prec <- matrix(c(2, 0.5, 0.5, 2), nrow = 2)
  R <- cholesky(prec)
  F_mat <- matrix(c(1, 0, 0, 1), nrow = 2)  # x >= 0, y >= 0
  g_vec <- c(0, 0)
  init <- c(1.5, 2.5)
  
  first <- harmonicHMC(
    nSample = 3,
    mean = mean_vec,
    choleskyFactor = R,
    constrainDirec = F_mat,
    constrainBound = g_vec,
    init = init,
    precFlg = TRUE,
    seed = 1
  )
  
  second <- harmonicHMC(
    nSample = 3,
    mean = mean_vec,
    choleskyFactor = R,
    constrainDirec = F_mat,
    constrainBound = g_vec,
    init = init,
    precFlg = TRUE,
    seed = 1
  )
  
  expect_equal(first, second)
})

test_that("harmonicHMC works with extraOutputs", {
  d <- 2
  mean_vec <- rep(0, d)
  prec <- diag(1, d)
  R <- cholesky(prec)
  F_mat <- diag(d)
  g_vec <- rep(0, d)
  init <- rep(1, d)
  
  result <- harmonicHMC(
    nSample = 3,
    mean = mean_vec,
    choleskyFactor = R,
    constrainDirec = F_mat,
    constrainBound = g_vec,
    init = init,
    precFlg = TRUE,
    seed = 456,
    extraOutputs = "numBounces"
  )
  
  expect_type(result, "list")
  expect_named(result, c("samples", "numBounces"))
  expect_equal(dim(result$samples), c(3, d))
  expect_length(result$numBounces, 3)
})

test_that("getHarmonicSample works with whitened coordinates", {
  # Create simple whitened constraints for testing
  whitened_constraints <- list(
    direc = matrix(c(1, 0, 0, 1), nrow = 2),
    direcRowNormSq = c(1, 1),
    bound = c(-0.5, -0.5)
  )
  
  whitened_pos <- c(0.1, 0.2)
  
  result <- getHarmonicSample(
    whitenedPosition = whitened_pos,
    whitenedConstraints = whitened_constraints,
    integrationTime = pi/4
  )
  
  expect_type(result, "list")
  expect_true("position" %in% names(result))
  expect_true("numBounces" %in% names(result))
  expect_length(result$position, 2)
  expect_type(result$numBounces, "integer")
})