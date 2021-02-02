#' Skewness (Copied from the deprecated R-package modes)
#' @export
#' @description Computes skewness.
#' @details Due to discontinuation from CRAN (Package ‘modes’ was removed from the CRAN repository; Archived on 2020-03-03 as check
#' problems were not corrected despite reminders), we have copied the required functions into UMtools to ease installation (and slightly modified).
#' The original documentation states the following: "This function calculates the skewness of a data with optional bias correction.
#' The skewness is a measure of the symmetry of a distrbution. A negative skewness means the data is left skewed or has a fat left tail.
#' The converse is true for a positive skew".
#'
#' @param x Data vector.
#' @param finite Should the finite sample size correction be applied? Defaults to TRUE.
#' @return Skewness (numeric)
#' @examples
#' data <- c(rnorm(15, 0, 1), rnorm(21, 5, 1))
#' hist(data)
#' skewness(data, TRUE)
skewness <- function(x, finite = TRUE)
{
  n = length(x)
  S = (1/n) * sum((x - mean(x))^3)/(((1/n) * sum((x - mean(x))^2))^1.5)
  if (finite == FALSE) {
    S = S
  }
  else {
    S = S * (sqrt(n * (n - 1)))/(n - 2)
  }
  return(S)
}
