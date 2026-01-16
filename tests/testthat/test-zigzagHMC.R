test_that("zigzagHMC matches TruncatedNormal reference for truncated Gaussian", {
  skip_if_not_installed("TruncatedNormal")
  library(TruncatedNormal)
  set.seed(123)
  
  # Define problem
  d <- 2
  meanVec <- c(0, 0)
  covMat <- matrix(c(1, 0.5, 0.5, 1), nrow = d)
  precMat <- solve(covMat)  # Convert to precision matrix
  lb <- c(-1, -1)
  ub <- c(Inf, Inf)
  
  nSamples <- 50000
  burnin <- 20000
  nRef <- 100000
  
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
  
  samples_zz <- zigzagHMC(
    nSample = nSamples,
    burnin = burnin,
    mean = meanVec,
    prec = precMat,
    lowerBounds = lb,
    upperBounds = ub,
    init = init,
    nutsFlg = FALSE,  
    seed = 123,
    diagnosticMode = FALSE
  )
  
  # Compare means
  ref_means <- colMeans(samples_ref)
  zz_means <- colMeans(samples_zz)
  
  # Print for debugging
  cat("\nReference means:", ref_means)
  cat("\nZigzag HMC means:", zz_means, "\n")
  
  tol <- 0.01
  expect_equal(zz_means, ref_means, tolerance = tol)
})

test_that("zigzagHMC produces consistent results with seeding", {
  d <- 2
  prec <- diag(c(1, 2))
  
  first <- zigzagHMC(
    nSample = 3, 
    mean = rep(0, d),
    prec = prec,
    lowerBounds = rep(0, d),
    upperBounds = rep(Inf, d),
    seed = 1
  )
  
  second <- zigzagHMC(
    nSample = 3,
    mean = rep(0, d),
    prec = prec,
    lowerBounds = rep(0, d),
    upperBounds = rep(Inf, d),
    seed = 1
  )
  
  expect_equal(first, second)
})

test_that("zigzagHMC works with different nutsFlg options", {
  d <- 2
  
  # Zigzag-HMC (nutsFlg = FALSE)
  samples_hmc <- zigzagHMC(
    nSample = 3,
    mean = rep(0, d),
    prec = diag(1, d),
    lowerBounds = rep(-Inf, d),
    upperBounds = rep(Inf, d),
    nutsFlg = FALSE,
    seed = 42
  )
  expect_equal(dim(samples_hmc), c(3, d))
  
  # Zigzag-NUTS (nutsFlg = TRUE)
  samples_nuts <- zigzagHMC(
    nSample = 3,
    mean = rep(0, d),
    prec = diag(1, d),
    lowerBounds = rep(-Inf, d),
    upperBounds = rep(Inf, d),
    nutsFlg = TRUE,
    seed = 42
  )
  expect_equal(dim(samples_nuts), c(3, d))
})

test_that("zigzagHMC diagnostic mode works", {
  d <- 2
  result <- zigzagHMC(
    nSample = 3,
    mean = rep(0, d),
    prec = diag(1, d),
    lowerBounds = rep(-Inf, d),
    upperBounds = rep(Inf, d),
    diagnosticMode = TRUE,
    seed = 456
  )
  
  expect_type(result, "list")
  expect_named(result, c("samples", "stepsize"))
  expect_equal(dim(result$samples), c(3, d))
})