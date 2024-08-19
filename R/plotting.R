utils::globalVariables(c("tau", "value", "aes","value", "segment"))

#' Plot method for decp_result
#' 
#' @param x An object of class 'decp_result'
#' @param ... Additional arguments passed to the plotting function
#' @export
plot.decp_result <- function(x, ...) {
  all_data <- data.frame(value = numeric(0), segment = factor())
  
  # Convert estimation_results into a data frame for plotting
  for (i in seq_along(x$estimation_results)) {
    vec_name <- paste("estimation_vecRW_", sprintf("%02d", i), sep = "")
    shifted_vector <- x$estimation_results[[vec_name]] + x$ordered_change_points[i]
    
    temp_df <- data.frame(value = shifted_vector,
                          segment = as.factor(paste("Segment", i)))
    all_data <- rbind(all_data, temp_df)
  }
  
  # Create the plot using ggplot2
  p <- ggplot2::ggplot(all_data, ggplot2::aes(x = value, fill = segment, color = segment)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::labs(x = "Time Series", y = "Density") +
    ggplot2::theme_minimal()
  
  return(p)
}


#' Plot method for mle_change_point_result
#' 
#' @param x An object of class 'mle_change_point_result'
#' @param ... Additional arguments passed to the plotting function
#' @export
plot.mle_change_point_result <- function(x, ...) {
  plot_data <- x$output_data
  best_tau <- x$best_tau
  best_value <- x$best_value
  
  ggplot2::ggplot(plot_data, ggplot2::aes(x = tau, y = value)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(ggplot2::aes(color = (tau == best_tau)), size = 2) +
    ggplot2::labs(title = "MLE Change Point Detection",
                  x = "Time Series",
                  y = "MLE value") +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_manual(values = c("black", "lightsalmon2"), guide = "none") +
    ggplot2::geom_vline(xintercept = best_tau, linetype = "dashed", color = "blue4") +
    ggplot2::annotate("text", x = best_tau, y = best_value, 
                      label = paste("Best Tau =", best_tau), vjust = -1.5)
}
