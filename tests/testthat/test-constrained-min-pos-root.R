test_that("first positive time is calculated correctly", {
  intercept <- -0.123
  slope <- 4.56
  expect_equal(firstPositiveTime(intercept, slope), - intercept / slope)
  
  slope <- -4.56
  expect_equal(firstPositiveTime(intercept, slope), Inf)
  
  intercept <- 1.23
  expect_equal(firstPositiveTime(intercept, slope), 0)
  
  slope <- -0.456
  expect_equal(firstPositiveTime(intercept, slope), 0)
})


test_that("constrained quadratic equation is solved correctly", {
  
  quadratic_func <- function (a, b, c, x) { a * x^2 + b * x + c }
  
  min_val <- -.123
  min_loc <- .456
  a <- .789
  
  b <- - 2 * a * min_loc
  c <- a * min_loc^2 + min_val
  
  constrainedRoot <- hdtg:::minimumPositiveRootWithConstraint(a, b, c, lowerBd = 0)
  expect_lt(constrainedRoot, min_loc)
  expect_equal(quadratic_func(a, b, c, constrainedRoot), 0)
  
  constrainedRoot <- hdtg:::minimumPositiveRootWithConstraint(a, b, c, lowerBd = 0.5)
  expect_gt(constrainedRoot, lowerBd)
  expect_equal(quadratic_func(a, b, c, constrainedRoot), 0)
  
  constrainedRoot <- hdtg:::minimumPositiveRootWithConstraint(a, b, c, lowerBd = 1)
  expect_equal(constrainedRoot, Inf)
  
  min_loc <- 0
  b <- - 2 * a * min_loc
  c <- a * min_loc^2 + min_val
  
  constrainedRoot <- hdtg:::minimumPositiveRootWithConstraint(a, b, c, lowerBd = 0)
  expect_equal(quadratic_func(a, b, c, constrainedRoot), 0)
  
  constrainedRoot <- hdtg:::minimumPositiveRootWithConstraint(a, b, c, lowerBd = 1)
  expect_equal(constrainedRoot, Inf)
})