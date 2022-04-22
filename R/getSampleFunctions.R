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
