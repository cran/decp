#' MLE Change Point Detection
#'
#' @description Maximum likelihood estimation change point detection.
#' @param input_data A numeric matrix of observations for multivariate time series
#' data where the dimension is not greater than the observations. Date columns should not be inputted.
#' @param verbose Logical value indicating whether to print messages during the function execution. Default is TRUE.
#' @return An object of class 'mle_change_point_result' containing the index of the change point estimate, its MLE value, and the MLE data.
#' @export
#' @importFrom matrixcalc is.positive.definite
#' @examples
#' # Example usage
#' data <- matrix(rnorm(1000), ncol = 10)
#' tau_range <- 30:(nrow(data) - 30)
#' result <- mle_change_point_detection(data)
#' print(result)
mle_change_point_detection <- function(input_data, verbose = TRUE) {
  input_data <- data.matrix(input_data)
  
  # Function to calculate the value of the equation
  calculate_value <- function(dat, tau, num_iterations) {
    data1 <- dat[1:tau, ]
    data2 <- dat[(tau + 1):num_iterations, ]
    Sigma1 <- cov(data1)
    Sigma2 <- cov(data2)
    det_Sigma1 <- det(Sigma1)
    det_Sigma2 <- det(Sigma2)
    inv_Sigma1 <- solve(Sigma1)
    inv_Sigma2 <- solve(Sigma2)
    
    sum1 <- sum(apply(data1, 1, function(Rk) as.numeric(t(Rk) %*% inv_Sigma1 %*% Rk)))
    sum2 <- sum(apply(data2, 1, function(Rk) as.numeric(t(Rk) %*% inv_Sigma2 %*% Rk)))
    
    # Calculate the equation value
    value <- -tau * (log(det_Sigma1) + sum1) - (num_iterations - tau) * (log(det_Sigma2) + sum2)
    
    return(value)
  }
  
  # Number of observations
  num_iterations <- nrow(input_data)
  
  # Range for tau
  num_columns <- ncol(input_data)
  tau_dim <- max(4 * num_columns, 30)
  tau_range <- tau_dim:(num_iterations - tau_dim)
  
  # Initialize a vector to store values for each tau
  values <- numeric(length(tau_range))
  for (i in seq_along(tau_range)) {
    tau <- tau_range[i]
    values[i] <- calculate_value(input_data, tau, num_iterations)
  }
  
  best_tau_index <- which.max(values)
  best_tau <- tau_range[best_tau_index]
  best_value <- values[best_tau_index]
  
  if (verbose) {
    message("The change point is at the observation: ", best_tau, " with MLE-value: ", best_value)
  }
  
  output_data <- data.frame(tau = tau_range, value = values)
  
  result <- list(best_tau = best_tau, best_value = best_value, output_data = output_data)
  class(result) <- "mle_change_point_result"
  return(result)
}

#' Print method for 'mle_change_point_result' class
#'
#' @param x An object of class 'mle_change_point_result'.
#' @param ... Additional arguments (not used).
#' @export
print.mle_change_point_result <- function(x, ...) {
  cat("MLE Change Point Detection Result:\n")
  cat("Best tau (change point):", x$best_tau, "\n")
  cat("Best value (MLE):", x$best_value, "\n")
}
