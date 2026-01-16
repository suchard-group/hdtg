test_that("markovianZigzag matches TruncatedNormal reference for truncated Gaussian", {
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
  
  # markovianZigzag samples
  samples_mz <- markovianZigzag(
    nSample = nSamples,
    burnin = burnin,
    mean = meanVec,
    prec = precMat,
    lowerBounds = lb,
    upperBounds = ub,
    init = init,
    seed = 123,
    diagnosticMode = FALSE
  )
  
  # Compare means
  ref_means <- colMeans(samples_ref)
  mz_means <- colMeans(samples_mz)
  
  # Print for debugging
  cat("\nReference means:", ref_means)
  cat("\nMarkovian Zigzag means:", mz_means, "\n")
  
  tol <- 0.05
  expect_equal(mz_means, ref_means, tolerance = tol)
})

test_that("markovianZigzag produces consistent results with seeding", {
  d <- 2
  prec <- diag(c(1, 2))
  
  first <- markovianZigzag(
    nSample = 3,
    mean = rep(0, d),
    prec = prec,
    lowerBounds = rep(0, d),
    upperBounds = rep(Inf, d),
    seed = 1
  )
  
  second <- markovianZigzag(
    nSample = 3,
    mean = rep(0, d),
    prec = prec,
    lowerBounds = rep(0, d),
    upperBounds = rep(Inf, d),
    seed = 1
  )
  
  expect_equal(first, second)
})

test_that("getMarkovianZigzagSample works", {
  # Create engine
  engine <- createEngine(
    dimension = 2,
    lowerBounds = c(-1, -1),
    upperBounds = c(1, 1),
    seed = 123,
    mean = c(0, 0),
    precision = diag(2)
  )
  
  position <- c(0.1, -0.2)
  
  # Test with provided travel time
  result <- getMarkovianZigzagSample(
    position = position,
    engine = engine,
    travelTime = 0.5
  )
  
  expect_type(result, "list")
  expect_named(result, c("position", "velocity"))
  expect_length(result$position, 2)
  expect_length(result$velocity, 2)
  
  # Test with custom velocity
  custom_velocity <- c(1, -1)
  result2 <- getMarkovianZigzagSample(
    position = position,
    velocity = custom_velocity,
    engine = engine,
    travelTime = 0.5
  )
  
  expect_named(result2, c("position", "velocity"))
})