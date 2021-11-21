#' Title
#'
#' @param get_prec_product function that returns the precision matrix - vector product
#' @param mean a d-dimensional mean vector
#' @param position a d-dimensional vector of the initial position
#' @param momentum a d-dimensional vector of the initial momentum
#' @param t time length to simulate the Markov process
#'
#' @return
#'
#' @examples
hzz <- function(get_prec_product,
                mean,
                position,
                constraits,
                momentum,
                t) {
  debug_flg = F
  ndim = length(position)
  position <- .get_initial_position(position, constraits)
  velocity <- sign(momentum)
  gradient <- get_prec_product(position - mean)
  action <- get_prec_product(velocity + mean)
  
  # list to store useful info
  dynamics <-
    list(
      position = position,
      velocity = velocity,
      action = action,
      gradient = gradient,
      momentum = momentum,
      column = NULL,
      event_time = NULL,
      event_idx = NULL,
      event_type = NULL
    )
  
  time_remaining <- t
  while (time_remaining > 0) {
    first_bounce <-
      .get_next_bounce(
        dynamics$position,
        dynamics$velocity,
        dynamics$action,
        dynamics$gradient,
        dynamics$momentum,
        constraits
      )
    event_time <- first_bounce$t_event
    event_idx <- first_bounce$index
    event_type <- first_bounce$event_type
    if (debug_flg) {
      print(paste("event time", paste(event_time, collapse = " ")))
      print(paste("event type", paste(event_type, collapse = " ")))
      print(paste("position", paste(dynamics$position, collapse = " ")))
      print(paste("momentun", paste(dynamics$momentum, collapse = " ")))
      print(paste("gradient", paste(dynamics$gradient, collapse = " ")))
      print(paste("action", paste(dynamics$action, collapse = " ")))
      cat("\n")
    }
    
    if (time_remaining < event_time) {
      # No event during remaining time
      
      # update position
      dynamics$position <- dynamics$position + time_remaining * dynamics$velocity
      # update momentum
      # TODO: this update of momentum is only necessary when using NUTS
      half_time_squared = time_remaining * time_remaining / 2
      dynamics$momentum <-
        dynamics$momentum - time_remaining * dynamics$gradient - half_time_squared * dynamics$action
      
      time_remaining <- 0
    } else {
      dynamics <-
        modifyList(
          dynamics,
          list(
            column = get_prec_product(event_idx),
            event_time = event_time,
            event_idx = event_idx,
            event_type = first_bounce$event_type
          )
        )
      # update dynamics
      dynamics <- .update_dynamics(dynamics)
      # reflect velocity element
      dynamics$velocity[event_idx] <- -dynamics$velocity[event_idx]
      time_remaining <- time_remaining - event_time
    }
    if (debug_flg) {
      print("after update")
      print(paste("position", paste(dynamics$position, collapse = " ")))
      print(paste("momentun", paste(dynamics$momentum, collapse = " ")))
      print(paste("gradient", paste(dynamics$gradient, collapse = " ")))
      print(paste("action", paste(dynamics$action, collapse = " ")))
      cat("**************************")
      cat("\n")
    }
  }
  
  return(dynamics$position)
}


.get_next_bounce <-
  function(position,
           velocity,
           action,
           gradient,
           momentum,
           constraits) {
    grad_time <- .get_grad_time(action, gradient, momentum)
    bound_time <- .get_bound_time(position, velocity, constraits)
    
    idx_grad <- which.min(grad_time)
    idx_bound <- which.min(bound_time)
    grad_flg <- grad_time[idx_grad] < bound_time[idx_bound]
    
    type <- ifelse(grad_flg, "gradient", "boundary")
    idx <- ifelse(grad_flg, idx_grad, idx_bound)
    t_event <-
      ifelse(grad_flg, grad_time[idx_grad], bound_time[idx_bound])
    
    return(list(
      t_event = t_event,
      index = idx,
      event_type = type
    ))
  }


#' Title
#'
#' @param dynamics
#'
#' @return
#'
#' @examples
.update_dynamics <- function(dynamics) {
  p <- dynamics$position
  v <- dynamics$velocity
  a <- dynamics$action
  g <- dynamics$gradient
  m <- dynamics$momentum
  c <- dynamics$column
  time <- dynamics$event_time
  index <- dynamics$event_idx
  
  half_time_squared = time ^ 2 / 2
  two_v1 = 2 * v[index]
  
  p <- p + time * v
  m <- m - time * g - half_time_squared * a
  g <- g + time * a
  a <- a - two_v1 * c
  
  if (dynamics$event_type == "boundary") {
    # Reflect against binary boundary
    m[index] <- -m[index]
    p[index] <- 0
  }
  
  if (dynamics$event_type == "gradient") {
    # Set 0 momentum to avoid numeric error
    m[index] <- 0
  }
  return(modifyList(dynamics, list(
    position = p,
    momentum = m,
    gradient = g,
    action = a
  )))
}
