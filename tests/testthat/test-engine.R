test_that("createEngine works and returns expected structure", {
  engine <- createEngine(
    dimension = 2,
    lowerBounds = c(-1, -1),
    upperBounds = c(1, 1),
    seed = 123,
    mean = c(0, 0),
    precision = diag(2)
  )
  
  expect_type(engine, "list")
  expect_named(engine, "engine")
  expect_true(is(engine$engine, "externalptr"))
  expect_true(length(engine) == 1)
})

test_that("createNutsEngine works", {
  nuts_engine <- createNutsEngine(
    dimension = 2,
    lowerBounds = c(-2, -2),
    upperBounds = c(2, 2),
    seed = 456,
    stepSize = 0.1,
    mean = c(0.5, -0.5),
    precision = matrix(c(2, 0.3, 0.3, 2), nrow = 2)
  )
  
  expect_type(nuts_engine, "list")
  expect_named(nuts_engine, "engine")
  expect_true(is(nuts_engine$engine, "externalptr"))
})

test_that("setMean updates engine without error", {
  engine <- createEngine(
    dimension = 2,
    lowerBounds = c(-1, -1),
    upperBounds = c(1, 1),
    seed = 123,
    mean = c(0, 0),
    precision = diag(2)
  )
  
  # Should work without error
  expect_error(setMean(engine, mean = c(1, 2)), NA)
  
  # Test with different mean values
  expect_error(setMean(engine, mean = c(-1, 1)), NA)
  expect_error(setMean(engine, mean = c(0.5, -0.5)), NA)
})

test_that("setPrecision updates engine without error", {
  engine <- createEngine(
    dimension = 2,
    lowerBounds = c(-1, -1),
    upperBounds = c(1, 1),
    seed = 123,
    mean = c(0, 0),
    precision = diag(2)
  )
  
  # Test diagonal precision
  expect_error(setPrecision(engine, precision = diag(c(2, 3))), NA)
  
  # Test correlated precision
  corr_prec <- matrix(c(2, 0.5, 0.5, 2), nrow = 2)
  expect_error(setPrecision(engine, precision = corr_prec), NA)
  
  # Test identity again
  expect_error(setPrecision(engine, precision = diag(2)), NA)
})

test_that("engine functions integrate correctly", {
  engine <- createEngine(
    dimension = 2,
    lowerBounds = c(-1, -1),
    upperBounds = c(1, 1),
    seed = 123,
    mean = c(0, 0),
    precision = diag(2)
  )
  
  # Get multiple samples to ensure engine works
  positions <- matrix(0, nrow = 5, ncol = 2)
  current_pos <- c(0.1, 0.2)
  
  for (i in 1:5) {
    current_pos <- getZigzagSample(
      position = current_pos,
      nutsFlg = FALSE,
      engine = engine,
      stepSize = 0.1
    )
    positions[i, ] <- current_pos
    
    # Update parameters every other iteration
    if (i %% 2 == 0) {
      setMean(engine, mean = runif(2, -0.5, 0.5))
      setPrecision(engine, precision = diag(runif(2, 0.5, 2)))
    }
  }
  
  # All positions should be within bounds
  expect_true(all(positions >= -1 & positions <= 1))
  expect_true(!any(is.na(positions)))
})

test_that("engine works with getMarkovianZigzagSample", {
  engine <- createEngine(
    dimension = 2,
    lowerBounds = c(-1, -1),
    upperBounds = c(1, 1),
    seed = 123,
    mean = c(0, 0),
    precision = diag(2)
  )
  
  result <- getMarkovianZigzagSample(
    position = c(0.1, -0.2),
    engine = engine,
    travelTime = 0.5
  )
  
  expect_type(result, "list")
  expect_named(result, c("position", "velocity"))
  expect_length(result$position, 2)
  expect_length(result$velocity, 2)
  expect_type(result$position, "double")
  expect_type(result$velocity, "double")
})