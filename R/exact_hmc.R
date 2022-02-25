library(Rcpp)
library(profvis)

sourceCpp(here("src", "EigenMatVec.cpp"))


# moved to RcppBounceTime
compute_bounce_time = function(position,
                               momentum,
                               constraint_direc,
                               constraint_bound) {
  # formula 2.23
  fa = constraint_direc %*% momentum
  fb = constraint_direc %*% position
  U = sqrt(fa ^ 2 + fb ^ 2)
  phi = atan2(-fa, fb)
  reachable_idxs = which(U > abs(constraint_bound))
  if (length(reachable_idxs) == 0) {  # no bounces will occur
    return(list("bounce_time" = Inf, "constraint_idx" = NA))
  }
  times = -phi[reachable_idxs] + acos(-constraint_bound[reachable_idxs] / U[reachable_idxs])
  min_time_idx = which.min(times)
  constraint_idx = reachable_idxs[min_time_idx]
  bounce_time = times[min_time_idx]
  return(list("bounce_time" = bounce_time, "constraint_idx" = constraint_idx))
}

# Moved to RcppHamiltonian()
simulate_hamiltonian = function(position, momentum, time) {
  # formula 2.10
  new_position = momentum * sin(time) + position * cos(time)
  new_momentum = momentum * cos(time) - position * sin(time)
  return(list("position" = new_position, "momentum" = new_momentum))
}

# Moved to RcppBounceMomentum
compute_bounce_momentum = function(position,
                                   momentum,
                                   constraint_direc,
                                   constraint_row_normsq,
                                   bounce_idx) {
  # formula 2.30

  alpha = (constraint_direc[bounce_idx, ] %*% momentum)[1,1] / # need scalar instead of 1x1 matrix
    constraint_row_normsq[bounce_idx]
  return(momentum - 2 * alpha * constraint_direc[bounce_idx, ])
}


generate_whitened_sample = function(initial_position,
                                    constraint_direc,
                                    constraint_row_normsq,
                                    bounds,
                                    total_time) {
  position = initial_position
  momentum = rnorm(ncol(constraint_direc))
  travelled_time = 0
  repeat {
    bounce = RcppBounceTime(position,
                            momentum,
                            constraint_direc,
                            bounds)
    if (bounce$bounce_time < total_time - travelled_time) {
      bounce_time = bounce$bounce_time
      hamiltonian = RcppHamiltonian(position, momentum, bounce_time)
      position = hamiltonian$position
      momentum = RcppBounceMomentum(
        position,
        hamiltonian$momentum,
        constraint_direc,
        constraint_row_normsq,
        bounce$constraint_idx
      )
      travelled_time = travelled_time + bounce_time
    } else {
      bounce_time = total_time - travelled_time
      hamiltonian = RcppHamiltonian(position, momentum, bounce_time)
      return(hamiltonian$position)
    }
  }
}


whiten_position = function(position,
                           constraint_direc,
                           constraint_bound,
                           cholesky,
                           mean,
                           paramet) {
  if (paramet == "prec") {
    return(cholesky %*% (position - mean))
  } else {
    return(backsolve(cholesky, position - mean, transpose=TRUE))
  }
}


unwhiten_position = function(position, cholesky, mean, paramet) {
  if (paramet == "prec") {
    return(backsolve(cholesky, position) + mean)
  } else {
    return(t(cholesky) %*% position + mean)
  }
}


whiten_constraints = function(constraint_direc,
                              constraint_bound,
                              cholesky,
                              mean,
                              paramet) {
  if (paramet == "prec") {
    direc = t(backsolve(cholesky, t(constraint_direc), transpose=TRUE))
    return(
      list(
        "direc" = direc,
        "direc_rownorm_sq" = rowSums(direc**2),
        "bound" = constraint_bound + constraint_direc %*% mean
      )
    )
  } else {
    direc = constraint_direc %*% t(cholesky)
    return(
      list(
        "direc" = direc,
        "direc_rownorm_sq" = rowSums(direc**2),
        "bound" = constraint_bound + constraint_direc %*% mean
      )
    )
  }
}


run_sampler_example = function(n,
                               initial_position,
                               constraint_direc,
                               constraint_bound,
                               cholesky,
                               mean,
                               paramet = c("prec", "cov"),
                               total_time = pi / 2,
                               seed = 1) {
  set.seed(seed)
  paramet = match.arg(paramet)
  whitened_constraints = whiten_constraints(constraint_direc,
                                            constraint_bound,
                                            cholesky,
                                            mean,
                                            paramet)
  sample = initial_position
  results = matrix(nrow = ncol(constraint_direc), ncol = n)
  for (i in 1:n) {
    sample = whiten_position(sample,
                             constraint_direc,
                             constraint_bound,
                             cholesky,
                             mean,
                             paramet)
    sample = generate_whitened_sample(
      sample,
      whitened_constraints$direc,
      whitened_constraints$direc_rownorm_sq,
      whitened_constraints$bound,
      total_time
    )
    sample = unwhiten_position(sample, cholesky, mean, paramet)
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
d = 100  # tested up to 1000, which generates 20 samples in 70s
A = matrix(runif(d^2)*2-1, ncol=d)
Sigma = t(A) %*% A
#Sigma = diag(d)
M = solve(Sigma)
mu = rep(0,d)
constraint_direc = diag(d)
constraint_bound = rep(0,d)

# check I didn't mess up dimensions
stopifnot(length(constraint_bound) == nrow(constraint_direc))

# run sampler in covariance mode
#profvis({
ptm = proc.time()
results = run_sampler_example(
  10000,
  rep(1, d),
  constraint_direc,
  constraint_bound,
  cholesky = chol(Sigma),
  mean = mu,
  paramet = "cov"
)
proc.time() - ptm
#})
rowMeans(results)
#var(t(results))

# run sampler in precision mode
ptm = proc.time()
results = run_sampler_example(
  10000,
  rep(1, d),
  constraint_direc,
  constraint_bound,
  cholesky = chol(M),
  mean = mu,
  paramet = "prec"
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
