#' @name plot_change_points
#' @title Plot Change Points
#' @description This function creates a density plot of change points from estimation results.
#' 
#' @param estimation_results A list of estimation results.
#' @param ordered_change_points A vector of ordered change points.
#' @return A ggplot object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_density labs
#' @importFrom utils globalVariables
#' @examples
#' # Example usage
#' estimation_results <- list(estimation_vecRW_01 = rnorm(100), estimation_vecRW_02 = rnorm(100))
#' ordered_change_points <- c(50, 150)
#' plot_change_points(estimation_results, ordered_change_points)
utils::globalVariables(c("Value", "Segment"))

plot_change_points <- function(estimation_results, ordered_change_points) {
  all_data <- data.frame(Value = numeric(0), Segment = factor())
  for (i in seq_along(estimation_results)) {
    vec_name <- paste("estimation_vecRW_", sprintf("%02d", i), sep = "")
    shifted_vector <- estimation_results[[vec_name]] + ordered_change_points[i]
    
    temp_df <- data.frame(Value = shifted_vector,
                          Segment = as.factor(paste("Segment", i)))
    all_data <- rbind(all_data, temp_df)
  }
  
  p <- ggplot(all_data, aes(x = Value, fill = Segment, color = Segment)) +
    geom_density(alpha = 0.5) +
    labs(x = "Time Series", y = "Density")
  
  return(p)
}

#' @name plot_mle_change_point
#' @title Plot MLE Change Point
#' @description This function creates a plot of MLE change points.
#' 
#' @param plot_data A data frame containing 'tau' and 'value' columns.
#' @param best_tau The best tau value for the change point.
#' @param best_value The maximum MLE value.
#' @return A ggplot object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_vline annotate
#' @importFrom utils globalVariables
#' @examples
#' # Example usage
#' plot_data <- data.frame(tau = 1:100, value = rnorm(100))
#' best_tau <- 50
#' best_value <- max(plot_data$value)
#' plot_mle_change_point(plot_data, best_tau, best_value)
utils::globalVariables(c("tau", "value"))

plot_mle_change_point <- function(plot_data, best_tau, best_value) {
  p <- ggplot(plot_data, aes(x = tau, y = value)) +
    geom_line() +
    geom_point(aes(color = (tau == best_tau)), size = 2) +
    geom_vline(xintercept = best_tau, linetype = "dashed") +
    annotate("text", x = best_tau, y = best_value,
             label = paste("Best Tau =", best_tau), vjust = -1.5)
  
  return(p)
}
