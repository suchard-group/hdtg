library(here)
library(profvis)
library(Rcpp)

sourceCpp(here("src", "ExactHMC.cpp"))


run_sampler_example = function(n,
                               initial_position,
                               constraint_direc,
                               constraint_bound,
                               cholesky,
                               mu,
                               precision = TRUE,
                               total_time = pi / 2,
                               seed = 1) {
  set.seed(seed)
  results = matrix(nrow = ncol(constraint_direc), ncol = n)
  whitened_constraints = WhitenConstraints(constraint_direc,
                                           constraint_bound,
                                           cholesky,
                                           mu,
                                           precision)
  sample = initial_position
  for (i in 1:n) {
    initial_momentum = rnorm(ncol(constraint_direc))
    sample = GenerateSample(
      sample,
      initial_momentum,
      whitened_constraints$direc,
      whitened_constraints$direc_rownorm_sq,
      whitened_constraints$bound,
      cholesky,
      mu,
      total_time,
      precision
    )
    results[, i] = sample
  }
  return(results)
}

# Example 1
d = 2 # dimension of independent multivariate normal

constraint_direc = matrix(c(1, 1, 1, 0, 0, 1), ncol = d, byrow = TRUE)  # works
constraint_bound = c(0, 0.5,-0.5)

# Example 2
#constraint_direc = matrix(c(1,1), ncol=d, byrow=TRUE)  # works
#constraint_bound = c(0)

# Example 3
#constraint_direc = matrix(c(1, 1, 1, 0), ncol = d, byrow = TRUE)  # works
#constraint_bound = c(0, 0.5)

#Example 4
# d = 1
# constraint_direc = matrix(c(1, 1), ncol = d, byrow = TRUE)  # works
# constraint_bound = c(0.1)


# # Example 5
# d = 4
# constraint_direc = matrix(c(1, 1, 1, 1,
#                             0, 1, 0, 0,
#                             0, 0, 1, 0,
#                             0, 0, 1, 1),
#                           ncol = d,
#                           byrow = TRUE)  # works
# constraint_bound = c(-1,-1,-0.5, 0)

# # Standard normal examples
# Sigma = diag(d)
# mu = rep(0,d)
# M = solve(Sigma)
# R = chol(M)



# Example 6: use non identity covariance and nonzero mean
d = 2
Sigma = matrix(c(10, 3, 3, 2), 2, 2)
mu = matrix(c(0.5, 1), ncol=1)
M = solve(Sigma)


# Example 7:
set.seed(1)
d = 100
A = matrix(runif(d^2)*2-1, ncol=d)
Sigma = t(A) %*% A
#Sigma = diag(d)
M = solve(Sigma)
mu = rep(0,d)
constraint_direc = diag(d)
constraint_direc = rbind(constraint_direc, c(1,1, rep(0,d-2)))
constraint_bound = runif(d+1, -0.2,0.2)

# check I didn't mess up dimensions
stopifnot(length(constraint_bound) == nrow(constraint_direc))

# run sampler in covariance mode
#profvis({
ptm = proc.time()
R = chol(Sigma)
results = run_sampler_example(
  10000,
  rep(1, d),
  constraint_direc,
  constraint_bound,
  cholesky = R,
  mu = mu,
  precision = FALSE
)
proc.time() - ptm
#})
rowMeans(results)

# verify constraints
dim(matrix(results[, which(colSums((constraint_direc %*% results) + constraint_bound < 0) == nrow(constraint_direc))], nrow = d))


#var(t(results))

# run sampler in precision mode
ptm = proc.time()
R = chol(M)
results = run_sampler_example(
  10000,
  rep(1, d),
  constraint_direc,
  constraint_bound,
  cholesky = R,
  mu = mu,
  precision = TRUE
)
proc.time() - ptm
rowMeans(results)
#var(t(results))


# # simulate naive way for verification
# library(MASS)
# ptm = proc.time()
# X = t(mvrnorm(n = 5000000, mu, Sigma))
# # only keep samples which satisfy all constraints
# X = matrix(X[, which(colSums((constraint_direc %*% X) + constraint_bound >= 0) == nrow(constraint_direc))], nrow = d)
# proc.time() - ptm
# ncol(X)  # number of valid samples
# rowMeans(X)
# var(t(X))
