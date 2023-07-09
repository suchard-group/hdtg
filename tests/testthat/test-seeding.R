test_that("results remain the same across runs", {
  d <- 4
  mode <- rep(0, d)
  sigma <- diag(0.5, d) + .5
  prec <- solve(sigma)
  
  firstSamples <- zigzagHMC(
    nSample = 5, mean = mode, prec = prec,
    lowerBounds = rep(0, d), 
    upperBounds = rep(Inf, d), 
    nutsFlg = TRUE,
    seed = 1
  )
  secondSamples <- zigzagHMC(
    nSample = 5, mean = mode, prec = prec,
    lowerBounds = rep(0, d), 
    upperBounds = rep(Inf, d), 
    nutsFlg = TRUE,
    seed = 1
  )
  expect_equal(firstSamples, secondSamples)
})
