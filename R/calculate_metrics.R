#' Calculate Metrics
#' 
#' @param lambda1 Eigenvalues of the first segment.
#' @param lambda2 Eigenvalues of the second segment.
#' @return A list of calculated metrics.
#' @keywords internal
#' @noRd
calculate_metrics <- function(lambda1, lambda2) {
  htta <- (lambda2 - lambda1) / lambda2
  zhta1 <- sqrt(sum((htta)^2))
  zhta2 <- sqrt(sum(((lambda2 - lambda1) / lambda1)^2))
  
  max_zhta <- max(zhta1, zhta2)
  
  x <- (lambda2 - lambda1) / lambda2
  y <- (lambda2 - lambda1) / lambda1
  
  x <- ifelse(x < 1, x, 1 - x)
  y <- ifelse(y > -1, y, 1 - y)
  
  term1 <- sum(log(1 - x) + x)
  term2 <- sum(log(1 + y) - y)
  
  return(list(htta = htta, zhta1 = zhta1, zhta2 = zhta2, max_zhta = max_zhta, term1 = term1, term2 = term2))
}
