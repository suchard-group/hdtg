# ' Generate one sample with ZZ-HMC
# '
# ' @param position a d-dimensional vector of the initial position
# ' @param momentum a d-dimensional vector of the initial momentum
# ' @param nutsFlg logical. `TRUE` for using No-U-Turn algorithm.
# ' @param engine ZZ-HMC engine.
# ' @param t time length to simulate the Markov process
# '
# ' @return one MCMC sample from the target MTN
getSample <- function(position,
                       momentum,
                       t,
                       nutsFlg,
                       engine) {
  if(is.null(momentum)){
    momentum <- drawMomentum(length(position))
  }
  
  if (nutsFlg) {
    res <- .oneNutsIteration(
      sexp = engine$engine,
      position = position,
      momentum = momentum
    )
    
  } else {
    res <- .oneIteration(
      sexp = engine$engine,
      position = position,
      momentum = momentum,
      time = t
    )
  }
  return(res$position)
}

drawMomentum <- function(dim){
  return((2 * (stats::runif(dim) > .5) - 1) * stats::rexp(dim, rate = 1))
}
