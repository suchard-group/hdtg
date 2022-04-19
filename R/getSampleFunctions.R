#' Title
#'
#' @param position a d-dimensional vector of the initial position
#' @param momentum a d-dimensional vector of the initial momentum
#' @param nutsFlg 
#' @param engine 
#' @param t time length to simulate the Markov process
#'
#' @return
#' @export
#' @examples
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
      momentum = momentum,
      stepsize = t
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

#' Title
#'
#' @param energyGrad
#' @param mean
#' @param position
#' @param constraits
#' @param momentum
#' @param t
#'
#' @return
#' @export
#'
#' @examples
getSampleR <- function(position,
                         momentum,
                         t,
                         energyGrad,
                         mean,
                         constraits) {
  debug_flg <- F
  
  velocity <- sign(momentum)
  gradient <- energyGrad(position - mean)
  action <- energyGrad(velocity)
  
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
      dynamics$position <-
        dynamics$position + time_remaining * dynamics$velocity
      # update momentum
      # TODO: this update of momentum is only necessary when using NUTS
      half_time_squared <- time_remaining * time_remaining / 2
      dynamics$momentum <-
        dynamics$momentum - time_remaining * dynamics$gradient - half_time_squared * dynamics$action
      
      time_remaining <- 0
    } else {
      dynamics$column <- energyGrad(event_idx)
      dynamics$event_time <- event_time
      dynamics$event_idx <- event_idx
      dynamics$event_type <- event_type
      # update dynamics
      dynamics <- .update_dynamics(dynamics)
      # reflect velocity element
      dynamics$velocity[event_idx] <-
        -dynamics$velocity[event_idx]
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


.update_dynamics <- function(dy) {
  half_time_squared <- dy$event_time ^ 2 / 2
  two_v1 = 2 * dy$velocity[dy$event_idx]
  
  dy$position <- dy$position + dy$event_time * dy$velocity
  dy$momentum <-
    dy$momentum - dy$event_time * dy$gradient - half_time_squared * dy$action
  dy$gradient <- dy$gradient + dy$event_time * dy$action
  dy$action <- dy$action - two_v1 * dy$column
  
  if (dy$event_type == "boundary") {
    # Reflect against binary boundary
    dy$momentum[dy$event_idx] <- -dy$momentum[dy$event_idx]
    dy$position[dy$event_idx] <- 0
  }
  
  if (dy$event_type == "gradient") {
    # Set 0 momentum to avoid numeric error
    dy$momentum[dy$event_idx] <- 0
  }
  return(dy)
}
