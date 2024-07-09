#' Simulate Estimation
#' 
#' @description The estimation of the detected change point.
#' @param lambda1 Eigenvalues of the first segment.
#' @param lambda2 Eigenvalues of the second segment.
#' @param term1 The negative drift term of the left hand side of the random walk.
#' @param term2 The negative drift term of the right hand side of the random walk.
#' @param num_iterations Determines the size of the two-sided random walk in the estimation process (each path).
#' If the jump size of the change point is small, num_iterations should be set to higher values to achieve accurate results.
#' For jump size >= 1, the default value is 100.
#' @param num_simulations Specifies the number of simulations to be conducted during the estimation process.
#' It is recommended to set num_simulations to a large value to ensure greater certainty and reliability of the results.
#' A higher number of simulations helps in capturing the variability and improves the accuracy of the estimation.
#' @return A numeric vector of the estimation results centered around zero.
#' The spike of the histogram is represents estimated change point, and it is expected to be at zero.
#' @export
#' @importFrom stats rnorm
#' @examples
#' # Example usage
#' lambda1 <- rnorm(10)
#' lambda2 <- rnorm(10)
#' term1 <- -1
#' term2 <- -2
#' result <- simulate_estimation(lambda1, lambda2, term1, term2, 
#'                                num_iterations = 100, num_simulations = 100)
simulate_estimation <- function(lambda1, lambda2, term1, term2,
                               num_simulations, num_iterations) {
  num_dimensions <- length(lambda1)
  estimation_vec1 <- numeric(num_simulations)
  
  for (k in 1:num_simulations) {
    left <- right <- numeric(num_iterations)
    
    for (t in 1:num_iterations ) {
      ra <- rnorm(num_dimensions)
      rb <- rnorm(num_dimensions)
      
      left[t] <- sum((ra^2 - 1) * ((lambda2 - lambda1) / lambda2)) + term1
      right[t] <- -sum((rb^2 - 1) * ((lambda2 - lambda1) / lambda1)) + term2
    }
    
    p1 <- 0.5 * cumsum(left)
    p2 <- 0.5 * cumsum(right)
    
    if (max(max(p1), max(p2), 0) == 0) {
      estimation_vec1[k] <- 0
    } else if (max(p2) < max(p1)) {
      estimation_vec1[k] <- -which.max(p1)
    } else {
      estimation_vec1[k] <- which.max(p2)
    }
  }
  
  return(estimation_vec1)
}
