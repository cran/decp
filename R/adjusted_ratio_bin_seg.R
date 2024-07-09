#' Adjusted Ratio Binary Segmentation
#' 
#' @description Adjusted ratio binary segmentation.
#' @param input_data A numeric matrix of observations for multivariate time series
#' data where the dimension is not greater than the observations. Date columns should not be inputted.
#' @param minseglen Minimum segment length for detecting change points.
#' @param alpha Level of significance for calculating the confidence intervals.
#' @return A list with change points and segments.
#' @export
#' @importFrom purrr map_dbl map map_if map_lgl pmap_dbl
#' @importFrom stats qnorm
#' @examples
#' # Example usage
#' data <- matrix(rnorm(1000), ncol = 10)
#' result <- adjusted_ratio_bin_seg(data, minseglen = 30, alpha = 0.05)

adjusted_ratio_bin_seg <- function(input_data, minseglen, alpha) {
  input_data <- data.matrix(input_data)
  recursive_bin_seg <- function(input_data, start, end, minseglen, alpha,
                                change_points, segments) {
    segment_length <- end - start + 1
    if (segment_length <= 2 * minseglen) {
      return(list(change_points = change_points, segments = segments))
    }
    
    test_stats <- matrix_dist_test_stat(input_data[start:end, , drop = FALSE], minseglen)
    threshold <- qnorm(1 - alpha / 2)
    
    max_test_stat <- max(test_stats, na.rm = TRUE)
    tau_star <- which.max(test_stats)
    
    tau_star_adjusted <- tau_star + start - 1
    
    if (max_test_stat > threshold) {
      change_points <- c(change_points, tau_star_adjusted)
      segments <- rbind(segments, data.frame(start = start,
                                             end = tau_star_adjusted),
                        data.frame(start = tau_star_adjusted + 1, end = end))
      
      result_left <- recursive_bin_seg(input_data, start, tau_star_adjusted,
                                       minseglen, alpha, change_points, segments)
      change_points <- result_left$change_points
      segments <- result_left$segments
      
      result_right <- recursive_bin_seg(input_data, tau_star_adjusted + 1,
                                        end, minseglen, alpha, change_points, segments)
      return(list(change_points = result_right$change_points, segments = result_right$segments))
    } else {
      return(list(change_points = change_points, segments = segments))
    }
  }
  
  initial_segments <- data.frame(start = 1, end = nrow(input_data))
  results <- recursive_bin_seg(input_data, 1, nrow(input_data), minseglen,
                               alpha, numeric(0), initial_segments)
  
  results$segments <- results$segments[-1, ]
  return(results)
}
