library(wordspace)


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
  if (length(reachable_idxs) == 0) {
    # no bounces will occur
    return(list("bounce_time" = NA, "constraint_idx" = NA))
  }
  times = -phi[reachable_idxs] + acos(-constraint_bound[reachable_idxs] / U[reachable_idxs])
  min_time_idx = which.min(times)
  constraint_idx = reachable_idxs[min_time_idx]
  bounce_time = times[min_time_idx]
  return(list("bounce_time" = bounce_time, "constraint_idx" = constraint_idx))
}


simulate_hamiltonian = function(position, momentum, time) {
  # formula 2.10
  new_position = momentum * sin(time) + position * cos(time)
  new_momentum = momentum * cos(time) - position * sin(time)
  return(list("position" = new_position, "momentum" = new_momentum))
}


compute_bounce_momentum = function(position,
                                   momentum,
                                   constraint_direc,
                                   constraint_row_normsq,
                                   bounce_idx) {
  # formula 2.30
  alpha = as.numeric(# need scalar instead of 1x1 matrix
    constraint_direc[bounce_idx,] %*% momentum / constraint_row_normsq[bounce_idx])
  return(momentum - 2 * alpha * constraint_direc[bounce_idx,])
}


generate_sample = function(initial_position,
                           constraint_direc,
                           constraint_row_normsq,
                           bounds,
                           total_time = pi / 2) {
  position = initial_position
  momentum = rnorm(ncol(constraint_direc))
  travelled_time = 0
  repeat {
    bounce = compute_bounce_time(position, 
                                 momentum, 
                                 constraint_direc, 
                                 bounds)
    if (!is.na(bounce$bounce_time) &&
        bounce$bounce_time < total_time - travelled_time) {
      bounce_time = bounce$bounce_time
      hamiltonian = simulate_hamiltonian(position, momentum, bounce_time)
      position = hamiltonian$position
      momentum = compute_bounce_momentum(
        position,
        hamiltonian$momentum,
        constraint_direc,
        constraint_row_normsq,
        bounce$constraint_idx
      )
      travelled_time = travelled_time + bounce_time
    } else {
      bounce_time = total_time - travelled_time
      hamiltonian = simulate_hamiltonian(position, momentum, bounce_time)
      return(hamiltonian$position)
    }
  }
}



run_sampler = function(n,
                       constraint_direc,
                       constraint_bound,
                       seed = 1) {
  set.seed(seed)
  constraint_row_normsq = wordspace::rowNorms(constraint_direc) ** 2
  sample = rep(1, ncol(constraint_direc))  # hardcoded initial position for testing
  results = matrix(nrow = ncol(constraint_direc), ncol = n)
  for (i in 1:n) {
    sample = generate_sample(sample,
                             constraint_direc,
                             constraint_row_normsq,
                             constraint_bound)
    results[, i] = sample
  }
  return(results)
}


d = 2 # dimension of independent multivariate normal
# Example 1
constraint_direc = matrix(c(1,1,1,0,0,1), ncol=d, byrow=TRUE)  # works
constraint_bound = c(0, 0.5, -0.5)

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

# check I didn't mess up dimensions
stopifnot(length(constraint_bound) == nrow(constraint_direc))



# run sampler
ptm = proc.time()
results = run_sampler(100000, constraint_direc, constraint_bound)
proc.time() - ptm
rowMeans(results)
apply(results, 1, var)


# simulate naive way for verification
set.seed(1)
n = 5000000
ptm = proc.time()
X = matrix(rnorm(n * d), nrow = d)
# only keep samples which satisfy all constraints
X = matrix(X[, which(colSums((constraint_direc %*% X) + constraint_bound >= 0) == nrow(constraint_direc))], nrow = d)
proc.time() - ptm
ncol(X)  # valid samples
rowMeans(X)
apply(X, 1, var)