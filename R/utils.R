#' Get an initial position for HZZ
#'
#' @param p0 a d-dimensional vector of the initial value. It must satisfy all constraints. If not specified a random initial value will be used
#' @param constraits a (list?) of constraints
#' @NoRd
#' @return a d-dimensional vector of the initial value
#'
#' @examples
.get_initial_position <- function(p0, constraits){
  stopifnot('random generated initial position not implemented yet' = !is.null(p0))
  #if (!is.null(p0)){
  return(p0)
  #} else {
    # TODO generate a random initial position
  #}
}


#' To find the minimal positive root of ax^2 + bx + c = 0 if it exists. Otherwise return Inf. 
#'
#' @param a vector of a
#' @param b vector of b
#' @param c vector of c
#'
#' @return vector of the minimal positive roots
#'
#' @examples
.min_pos_root <- function(a, b, c) {
  
  sign_a <- sign(a)
  b <- b * sign_a
  c <- c * sign_a
  a <- abs(a)
  
  disc <- b ^ 2 - 4 * a * c
  
  root <- rep(Inf, length(a))
  has_root <- (disc > 0)
  
  sqrtdisc <- disc
  sqrtdisc[has_root] <- sqrt(sqrtdisc[has_root])

  root[has_root] = (-b[has_root] - sqrtdisc[has_root]) / (2 * a[has_root])
  
  idx_1 <- root <= 0.0
  root[idx_1] <- (-b[idx_1] + sqrtdisc[idx_1]) / (2 * a[idx_1])
  idx_2 <- root <= 0.0
  root[idx_2] <- Inf

  return(root)
}


#' Title
#'
#' @param action 
#' @param gradient 
#' @param momentum 
#'
#' @return
#'
#' @examples
.get_grad_time <- function(action, gradient, momentum){
  return(.min_pos_root(action / 2, gradient, -momentum))
}


#' Title
#'
#' @param position 
#' @param velocity 
#'
#' @return
#'
#' @examples
.get_bound_time <- function(position, velocity, constraits){
  times <- rep(Inf, length(position))
  towards_boundary <- constraits * velocity < 0
  times[towards_boundary] <- abs(position / velocity)[towards_boundary]
  return(times)
}

