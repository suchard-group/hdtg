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