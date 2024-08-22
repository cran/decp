#' Detect and Estimate Change Points
#' 
#' @description Detect and estimate change points.
#' @param input_data A numeric matrix of observations for multivariate time series
#' data where the dimension is not greater than the observations. Date columns should not be inputted.
#' @param alpha Level of significance for calculating the confidence intervals
#' @param num_iterations Determines the size of the two-sided random walk in the estimation process (each path).
#' If the jump size of the change point is small, num_iterations should be set to higher values to achieve accurate results.
#' For jump size >= 1, the default value is 100.
#' @param num_simulations Specifies the number of simulations to be conducted during the estimation process.
#' It is recommended to set num_simulations to a large value to ensure greater certainty and reliability of the results.
#' A higher number of simulations helps in capturing the variability and improves the accuracy of the estimation.
#' @param verbose Logical value indicating whether to print messages during the function execution. Default is TRUE.
#' @return An object of class 'decp_result' containing the ordered change points, the summary of the jump sizes for each pair of segments,
#' the Confidence Interval (C.I.) of each detected change point, the maximum zhta, confidence interval level, and warnings in case that the C.I. of two adjacent change points overlap.
#' @export
#' @importFrom purrr map_dbl map map_if map_lgl pmap_dbl
#' @importFrom stats cov qnorm na.omit quantile rnorm
#' @examples
#' # Example usage
#' data_part1 <- matrix(rnorm(1500, mean = 0, sd = 1), ncol = 5)
#' data_part2 <- matrix(rnorm(1500, mean = 3, sd = 1), ncol = 5)
#' data <- rbind(data_part1, data_part2)
#' result <- decp(data, alpha = 0.05, num_simulations = 100, num_iterations = 50)
#' print(result)
decp <- function(input_data, alpha = 0.05, num_simulations = 10000, num_iterations = 100, verbose = TRUE) {
  # Ensure data is in the correct format (numeric matrix)
  X <- data.matrix(input_data)
  
  # Apply R-K change point detection
  n <- nrow(X)
  p <- ncol(X)
  minseglen <- max(4 * p, 30)
  
  result <- matrix_dist_test_stat(X, minseglen - p)
  
  messages <- character()
  
  if (max(result, na.rm = TRUE) > qnorm(1 - alpha / n)) {
    detected_change_point <- which.max(abs(result))
  }
  
  results <- adjusted_ratio_bin_seg(X, minseglen = minseglen - p, alpha = alpha)
  ordered_change_points <- sort(results$change_points)
  
  message_text <- paste("The detected change points are:", paste(ordered_change_points, collapse = ", "))
  if (verbose) message(message_text)
  messages <- c(messages, message_text)
  
  if (length(ordered_change_points) < 2) {
    stop("Function Stopping: You need at least two change points to continue.")
  }
  
  segments <- list()
  start_point <- 1
  for (i in 1:length(ordered_change_points)) {
    segment_name <- paste0("segm_", sprintf("%02d", i))
    segments[[segment_name]] <- X[start_point:ordered_change_points[i], ]
    start_point <- ordered_change_points[i] + 1
  }
  segment_name <- paste0("segm_", sprintf("%02d", length(ordered_change_points) + 1))
  segments[[segment_name]] <- X[start_point:nrow(X), ]
  
  for (i in 1:length(segments)) {
    segment_name <- paste0("segm_", sprintf("%02d", i))
    assign(segment_name, as.data.frame(segments[[segment_name]]))
  }
  
  num_columns <- dim(X)[2]
  mu <- numeric(num_columns)
  
  for (i in 1:length(segments)) {
    segment_name <- paste0("segm_", sprintf("%02d", i))
    current_segment <- get(segment_name)
    if ("Date" %in% colnames(current_segment)) {
      current_segment <- current_segment[, !colnames(current_segment) %in% "Date", drop = FALSE]
    }
    if (!all(sapply(current_segment, is.numeric))) {
      next
    }
    
    covariance_matrix_name <- paste0("covariance_matrix_", segment_name)
    assign(covariance_matrix_name, cov(as.matrix(current_segment)))
    
    variance_name <- paste0("variance_", segment_name)
    assign(variance_name, diag(get(covariance_matrix_name)))
    
    lambda_name <- paste0("lambda_", sprintf("%02d", i))
    assign(lambda_name, eigen(get(covariance_matrix_name))$values)
  }
  
  all_metrics_results <- list()
  summary_max_zhta <- numeric(length(segments) - 1)
  
  for (i in 1:(length(segments) - 1)) {
    j <- i + 1
    lambda1_name <- paste0("lambda_", sprintf("%02d", i))
    lambda2_name <- paste0("lambda_", sprintf("%02d", j))
    
    lambda1 <- get(lambda1_name)
    lambda2 <- get(lambda2_name)
    
    metrics_results <- calculate_metrics(lambda1, lambda2)
    
    result_key <- paste("Segment", i, "vs Segment", j)
    all_metrics_results[[result_key]] <- metrics_results
    summary_max_zhta[i] <- metrics_results$max_zhta
  }
  
  if (verbose) {
    message("\nSummary:")
    for (k in 1:length(summary_max_zhta)) {
      message(paste("Pair (", k, ",", k + 1, ") - max zhta = ", summary_max_zhta[k]))
    }
  }
  
  ordered_change_points <- sort(unique(results$change_points))
  estimation_results <- list()
  
  for (i in 1:length(ordered_change_points)) {
    if (i == 1) {
      start_index <- 1
    } else {
      start_index <- ordered_change_points[i - 1] + 1
    }
    
    end_index <- ordered_change_points[i]
    segment_data <- X[start_index:end_index, ]
    
    if (i < length(ordered_change_points)) {
      next_segment_data <- X[(ordered_change_points[i] + 1):ordered_change_points[i + 1], ]
    } else {
      next_segment_data <- X[(ordered_change_points[i] + 1):nrow(X), ]
    }
    
    lambda1 <- eigen(cov(as.matrix(segment_data)))$values
    lambda2 <- eigen(cov(as.matrix(next_segment_data)))$values
    
    metrics_results <- calculate_metrics(lambda1, lambda2)
    
    result_key <- paste("Segment", i, "vs Segment", i + 1)
    all_metrics_results[[result_key]] <- metrics_results
    summary_max_zhta[i] <- metrics_results$max_zhta
    
    vec_name <- paste("estimation_vecRW_", sprintf("%02d", i), sep = "")
    estimation_vec1 <- simulate_estimation(lambda1, lambda2,
                                           metrics_results$term1, metrics_results$term2,
                                           num_iterations = num_iterations,
                                           num_simulations = num_simulations)
    estimation_results[[vec_name]] <- estimation_vec1
  }
  
  lower_percentile <- alpha / 2 * 100
  upper_percentile <- (1 - alpha / 2) * 100
  confidence_level <- (1 - alpha) * 100
  
  shifted_intervals <- list()
  warnings <- character()
  
  for (i in seq_along(estimation_results)) {
    vec_name <- paste("estimation_vecRW_", sprintf("%02d", i), sep = "")
    current_vector <- estimation_results[[vec_name]]
    cleaned_vector <- na.omit(current_vector)
    sorted_vector <- sort(cleaned_vector)
    ci <- quantile(sorted_vector, probs = c(lower_percentile / 100, upper_percentile / 100))
    
    shifted_lower_ci <- ci[1] + ordered_change_points[i]
    shifted_upper_ci <- ci[2] + ordered_change_points[i]
    shifted_intervals[[i]] <- c(shifted_lower_ci, shifted_upper_ci)
    
    if (verbose) {
      message(sprintf("For %s the %.0f%% C.I., of the change point %.0f, is [%.0f, %.0f]",
                      vec_name, confidence_level, ordered_change_points[i], shifted_lower_ci, shifted_upper_ci))
    }
  }
  
  for (i in 1:(length(shifted_intervals) - 1)) {
    current_interval <- shifted_intervals[[i]]
    next_interval <- shifted_intervals[[i + 1]]
    if (current_interval[2] >= next_interval[1]) {
      warning_text <- sprintf("Warning! Check the segment %d and %d.", i, i + 1)
      warnings <- c(warnings, warning_text)
      if (verbose) message(warning_text)
    }
  }
  
  if (length(warnings) == 0 && verbose) {
    message("No warnings.")
  }
  
  result_list <- list(
    estimation_results = estimation_results, 
    ordered_change_points = ordered_change_points,
    summary_max_zhta = summary_max_zhta,
    alpha = alpha,
    warnings = warnings
  )
  
  class(result_list) <- "decp_result"
  return(result_list)
}

#' Print method for decp_result
#' @param x An object of class 'decp_result'
#' @param ... Additional arguments (not used)
#' @export
print.decp_result <- function(x, ...) {
  cat("Detected Change Points:\n")
  cat("Number of change points:", length(x$ordered_change_points), "\n")
  cat("Change points:", paste(x$ordered_change_points, collapse = ", "), "\n")
}

#' Summary method for decp_result
#' @param object An object of class 'decp_result'
#' @param ... Additional arguments (not used)
#' @export
summary.decp_result <- function(object, ...) {
  cat("Summary of Change Point Detection:\n")
  cat("Number of change points detected:", length(object$ordered_change_points), "\n\n")
  
  cat("Change Points:\n")
  cat(paste(object$ordered_change_points, collapse = ", "), "\n\n")
  
  cat("Jump Sizes (max zhta for each segment pair):\n")
  for (k in 1:length(object$summary_max_zhta)) {
    cat(sprintf("Pair (%d, %d) - max zhta = %.4f\n", k, k + 1, object$summary_max_zhta[k]))
  }
  
  cat("\nConfidence Intervals:\n")
  
  alpha <- object$alpha  # Retrieve alpha from the result object
  lower_percentile <- alpha / 2 * 100
  upper_percentile <- (1 - alpha / 2) * 100
  confidence_level <- (1 - alpha) * 100
  
  for (i in seq_along(object$estimation_results)) {
    vec_name <- paste("estimation_vecRW_", sprintf("%02d", i), sep = "")
    current_vector <- object$estimation_results[[vec_name]]
    cleaned_vector <- na.omit(current_vector)
    sorted_vector <- sort(cleaned_vector)
    ci <- quantile(sorted_vector, probs = c(lower_percentile / 100, upper_percentile / 100))
    
    shifted_lower_ci <- ci[1] + object$ordered_change_points[i]
    shifted_upper_ci <- ci[2] + object$ordered_change_points[i]
    
    cat(sprintf("For %s the %.0f%% C.I., of the change point %.0f, is [%.0f, %.0f]\n",
                vec_name, confidence_level, object$ordered_change_points[i], shifted_lower_ci, shifted_upper_ci))
  }
  
  if (length(object$warnings) > 0) {
    cat("\nWarnings:\n")
    cat(paste(object$warnings, collapse = "\n"), "\n")
  } else {
    cat("\nNo warnings.\n")
  }
}
