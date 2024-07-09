#' Estimate Covariance Distance
#' 
#' @param data A matrix of data.
#' @param epsilon A distance parameter.
#' @param f A function to apply to the eigenvalues.
#' @return A vector of test statistics.
#' @keywords internal
#' @noRd
#' @importFrom purrr map accumulate map_dbl map_if map_lgl pmap_dbl
covariance_distance_estimator <- function(data, epsilon, f) {
  products <- map(as.data.frame(t(data)), ~ .x %*% t(.x))
  forward_cumsum <- accumulate(products, ~ .x + .y)
  backward_cumsum <- map(forward_cumsum, ~ -.x + forward_cumsum[[nrow(data)]])
  
  function(tau) {
    sigma1 <- (1 / tau) * forward_cumsum[[tau]] + epsilon * diag(ncol(data))
    sigma2 <- (1 / (nrow(data) - tau)) * backward_cumsum[[tau]] + epsilon * diag(ncol(data))
    
    output <- tryCatch({
      cov_dist(sigma1, sigma2, f)
    }, error = function(error_message) {
      return(NA)
    })
    
    return(output)
  }
}

#' Covariance Distance Calculation
#' 
#' @param sigma1 The first covariance matrix.
#' @param sigma2 The second covariance matrix.
#' @param f A function to apply to the eigenvalues.
#' @return The calculated covariance distance.
#' @keywords internal
#' @noRd
#' @importFrom rlang exec
cov_dist <- function(sigma1, sigma2, f) {
  A <- geigen::geigen(sigma2, sigma1, symmetric = TRUE)$values
  return(rlang::exec(f, (A - 1)^2) + rlang::exec(f, (1 / A - 1)^2))
}

#' Estimate Covariance Distance
#' 
#' @param data A matrix of data.
#' @param minseglen Minimum segment length.
#' @return A vector of test statistics.
#' @keywords internal
#' @noRd
#' @importFrom purrr map_dbl map map_if map_lgl pmap_dbl
#' @importFrom magrittr %>%
#' @importFrom rlang exec
matrix_dist_test_stat <- function(data, minseglen) {
  p <- ncol(data)
  n <- nrow(data)
  t <- seq(p + minseglen + 1, n - p - minseglen)
  estimate <- covariance_distance_estimator(data, 0, mean)
  test_stat <- map_dbl(t, estimate)
  dimension_over_length <- map(t, ~ c(p / .x, p / (n - .x)))
  trace <- map(dimension_over_length, ~ rlang::exec(calculate_expected_trace, !!!.x))
  values <- map_lgl(trace, ~ length(.x[[1]]) == 0)
  bias <- map_dbl(dimension_over_length, ~ rlang::exec(asymptotic_bias, !!!.x))
  variance <- map_dbl(dimension_over_length, ~ rlang::exec(asymptotic_variance, !!!.x))
  trace <- map_if(trace, values, ~ NA) %>% map_dbl(~ .x[[1]])
  bias <- map_if(bias, values, ~ NA)
  variance <- map_if(variance, values, ~ NA)
  test_stat <- pmap_dbl(list(test_stat, trace, bias, variance), ~(p * (..1 - ..2) - ..3) / sqrt(..4))
  return(c(rep(NA, p + minseglen), abs(test_stat), rep(NA, p + minseglen)))
}

#' @keywords internal
#' @noRd
function_prod <- function(f1, f2){
  return(function(x){ return(f1(x) * f2(x))})
}

#' @importFrom purrr safely
#' @importFrom stats integrate
#' @keywords internal
#' @noRd
calculate_expected_trace <- function(y1, y2){
  asymptotic_pdf <- fisher_esd(y1, y2)
  integrand <- function_prod(function(x){(1 - x)^2 + (1 - (1 / x))^2}, asymptotic_pdf)
  asymptotic_supports <- construct_fisher_support(y1, y2)
  safe_integral <- safely(integrate)
  integral <- exec(safe_integral, integrand, !!!asymptotic_supports)[[1]]
  return(integral)
}

#' @keywords internal
#' @noRd
asymptotic_bias <- function(y1, y2){
  h <- sqrt(y1 + y2 - y1 * y2)
  K_1 <- 2 * h * (1 + h^2) / (1 - y2)^4 - 2 * h / (1 - y2)^2
  J_1 <- 2 * h * (1 + h^2) / (1 - y1)^4 - 2 * h / (1 - y1)^2
  
  return(2 * (h^2 - y2^2) / (1 - y2)^4 + 2 * K_1 * y2 / h + 2 * (h^2 - y1^2) / (1 - y1)^4 + 2 * J_1 * y1 / h)
}

#' @keywords internal
#' @noRd
asymptotic_variance <- function(y1, y2){
  h <- sqrt(y1 + y2 - y1 * y2)
  
  K_21 <- 2 * h * (1 + h^2) / (1 - y2)^4 - 2 * h / (1 - y2)^2
  K_22 <- 2 * h * (1 + h^2) / (1 - y1)^4 - 2 * h / (1 - y1)^2
  K_31 <- h^2 / (1 - y2)^4
  K_32 <- h^2 / (1 - y1)^4
  J_1 <- -2 * (1 - y2)^2
  J_2 <- (1 - y2)^4
  
  var_x <- K_21^2 + 2 * K_31^2
  var_y <- K_22^2 + 2 * K_32^2
  
  cov_xy <- J_1 * K_21 / h + 
    J_1 * K_21 / (h * (h^2 - 1)) + 
    (-J_1 * K_31 * (h^2 + 1) / h^2) + 
    (-J_1 * K_31 / (h^2 * (h^2 - 1))) +
    J_2 * K_21 * 2 * h / (h^2 - 1)^3 +
    J_2 * K_31 / h^2 +
    J_2 * K_31 * ((1 - 3 * h^2) / (h^2 * (h^2 - 1)^3))
  
  return(3 * (var_x + var_y + 2 * cov_xy))
}

#' @keywords internal
#' @noRd
fisher_esd <- function(y1, y2) {
  c2 <- y2
  c1 <- y1
  h <- sqrt(c1 + c2 - c1 * c2)
  a <- ((1 - h)^2) / (1 - c2)^2
  b <- ((1 + h)^2) / (1 - c2)^2
  
  function(x) {
    result <- numeric(length(x))
    result[x > a & x < b] <- ((1 - c2) * sqrt((b - x) * (x - a))) / (2 * pi * x * (c1 + x * c2))
    return(result)
  }
}

#' @keywords internal
#' @noRd
construct_fisher_support <- function(y1, y2){
  c2 <- y2
  c1 <- y1
  h <- sqrt(c1 + c2 - c1 * c2)
  a <- ((1 - h)^2) / (1 - c2)^2
  b <- ((1 + h)^2) / (1 - c2)^2
  return(c(a, b))
}
