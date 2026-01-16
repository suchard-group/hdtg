test_that("drawLaplaceMomentum works", {
  # Test with different dimensions
  momentum3 <- drawLaplaceMomentum(3)
  momentum5 <- drawLaplaceMomentum(5)
  momentum1 <- drawLaplaceMomentum(1)
  
  expect_length(momentum3, 3)
  expect_length(momentum5, 5)
  expect_length(momentum1, 1)
  expect_type(momentum3, "double")
})

test_that("getInitialPosition works with different bounds", {
  # Test 1: All finite bounds
  pos1 <- getInitialPosition(
    mean = c(0, 0),
    lowerBounds = c(-1, -2),
    upperBounds = c(1, 2)
  )
  expect_length(pos1, 2)
  
  # Test 2: Mixed bounds
  pos2 <- getInitialPosition(
    mean = c(0, 0, 0),
    lowerBounds = c(-Inf, 0, -1),
    upperBounds = c(Inf, 5, Inf)
  )
  expect_length(pos2, 3)
  
  # Test 3: All infinite bounds (should return mean)
  pos3 <- getInitialPosition(
    mean = c(1, 2, 3),
    lowerBounds = c(-Inf, -Inf, -Inf),
    upperBounds = c(Inf, Inf, Inf)
  )
  expect_equal(pos3, c(1, 2, 3))
})

test_that("cholesky works and matches R's chol", {
  # Create a positive definite matrix
  set.seed(123)
  B <- matrix(rnorm(16), 4, 4)
  B <- t(B) %*% B  # Make symmetric positive definite
  
  # Use your cholesky function
  U <- cholesky(B)
  
  # Use R's built-in chol (returns upper triangular)
  U_R <- chol(B)
  
  expect_true(is.matrix(U))
  expect_equal(dim(U), c(4, 4))
  
  # Check it's upper triangular
  expect_true(all(U[lower.tri(U)] == 0))
  
  # Both should be valid Cholesky decompositions
  expect_equal(t(U) %*% U, B, tolerance = 1e-10)
  expect_equal(t(U_R) %*% U_R, B, tolerance = 1e-10)
  
  # Check that U and U_R produce same quadratic form
  test_vec <- rnorm(4)
  expect_equal(
    sum((t(U) %*% test_vec)^2),
    sum((t(U_R) %*% test_vec)^2),
    tolerance = 1e-10
  )
})