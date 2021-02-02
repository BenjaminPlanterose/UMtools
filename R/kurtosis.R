#' Excess kurtosis (Copied from the deprecated R-package modes)
#' @export
#' @description Computes excess kurtosis.
#' @details Due to discontinuation from CRAN (Package ‘modes’ was removed from the CRAN repository; Archived on 2020-03-03 as check
#' problems were not corrected despite reminders), we have copied the required functions into UMtools to ease installation (and slightly modified).
#' The original documentation states the following: "This function calculates the excess kurtosis of a data vector with optional
#' bias correction. Kurtosis is a meaure of the peakedness or how heavy the tails of a distribution are–this dual interpretation
#' is a result of the obvious inverse relationship between fat tails and high peaks. Excess kurtosis is simply "kurtosis-3."
#' This is a correction that is often done to allow for comparision to the normal distribution–which has a kurtosis of 3 and excess
#'  kurtosis of 0. A kurtosis greater than 0 means that the distribution is leptokurtic and so has a high peak with skinny tails.
#'  Conversely, a kurtosis less than 0 means that the distribution is platykurtic and so has a low peak and heavy tails.
#'  This interpretation is slightly more complicated once the distribution is not unimodal and/or non-zero skewness. Comparing to
#'  Gaussian (normal) moments is more acceptable in these cases".
#'
#' @param x Data vector.
#' @param finite Should the finite sample correction (bias correction) be used? Defaults to TRUE.
#' @return Kurtosis (numeric)
#' @examples
#' data <- c(rnorm(15, 0, 1), rnorm(21, 5, 1))
#' hist(data)
#' kurtosis(data, TRUE)
kurtosis <- function(x, finite)
{
  n = length(x)
  K = (1/n) * sum((x - mean(x))^4)/(((1/n) * sum((x - mean(x))^2))^2) -
    3
  if (finite == FALSE) {
    K = K
  }
  else {
    K = ((n - 1) * ((n + 1) * K - 3 * (n - 1))/((n - 2) * (n - 3))) + 3
  }
  return(K)
}
