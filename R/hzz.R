#' Title
#'
#' @param get_prec_product function that returns the precision matrix - vector product
#' @param mean a d-dimensional mean vector
#' @param position a d-dimensional vector of the initial position
#' @param momentum a d-dimensional vector of the initial momentum
#' @param t time length to simulate the Markov process
#'
#' @return
#' @export
#'
#' @examples
hzz <- function(get_prec_product, mean, position, momentum, t){
  
  ndim = length(position)
  
  position <- get_initial_position(position, constraits)
  momentum <-
    (2 * (runif(ndim) > .5) - 1) * rexp(ndim, rate = 1)
  velocity <- sign(momentum)
  gradient <- get_prec_product(position - mean)
  action <- get_prec_product(velocity + mean)
  
  time_remaining <- t
  #bounce_state <- initialize_state(t)

  while (time_remaining > 0) {
    first_bounce <- get_next_bounce(position, velocity, action, gradient, momentum)
    #bounce_state <- do_bounce(bounce_state, first_bounce, position, velocity, action, gradient, momentum)
    event_time <- first_bounce$t_event
    
    final_bounce_state
    if (time_remaining < eventTime) { # No event during remaining time
      
      updatePositionAndMomentum(position, velocity, action, gradient, momentum, remainingTime);
      # update position
      position <- position + time_remaining * velocity
      # update momentum
      # TODO: this update of momentum is only necessary when using NUTS
      half_time_squared = time_remaining * time_remaining / 2
      momentum <- momentum - remaining_time * gradient - half_time_squared * action

      time_remaining = 0
    } else {
      final Type eventType = firstBounce.type;
      final int eventIndex = firstBounce.index[0];
      
      WrappedVector column = getPrecisionColumn(eventIndex);
      
      updateDynamics(position, velocity, action, gradient, momentum, column, eventTime, firstBounce.index,
                     eventType);
      
      reflectVelocity(velocity, firstBounce.index);
      
      finalBounceState = new BounceState(eventType, eventIndex, remainingTime - eventTime);
      
      recordEvents(eventType);
    }
    
    if (TIMING) {
      timer.stopTimer("doBounce");
    }
    
    return finalBounceState;
  }

  return(position)
}


#' Title
#'
#' @param position 
#' @param velocity 
#' @param action 
#' @param gradient 
#' @param momentum 
#'
#' @return
#' @export
#'
#' @examples
get_next_bounce <- function(position, velocity, action, gradient, momentum){
  
  grad_time <- get_grad_time(action, gradient, momentum)
  bound_time <- get_bound_time(position, velocity)
  
  idx_grad <- which.min(grad_time)
  idx_bound <- which.min(bound_time)
  grad_flg <- grad_time[idx_grad] < bound_time[idx_bound]
  
  type <- ifelse(grad_flg, "gradient", "boundary")
  idx <- ifelse(grad_flg, idx_grad, idx_bound)
  t_event <- ifelse(grad_flg, grad_time[idx_grad], bound_time[idx_bound])
  
  return(list(t_event = t_event, index = idx, event_type = type))
}

