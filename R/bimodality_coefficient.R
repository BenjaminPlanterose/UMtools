#' Bimodality coefficient (Copied from the deprecated R-package modes)
#' @export
#' @description Computes bimodality coefficient.
#' @details Due to discontinuation from CRAN (Package ‘modes’ was removed from the CRAN repository; Archived on 2020-03-03 as check
#' problems were not corrected despite reminders), we have copied the required functions into UMtools to ease installation (and slightly modified).
#' The original documentation states the following: "This function calculates the Bimodality Coefficient of a data vector with the option
#' for a finite sample (bias) correction. This bias correction is important to correct for the (well-documented) finite sample bias.
#' The bimodality coefficient has a range of zero to one (that is: [0,1]) where a value greater than "5/9" suggests bimodality.
#' The maximum value of one ("1") can only be reached when the distribution is composed of two point masses".
#' @param x Data vector.
#' @param finite Should the finite sample size correction be applied to the skewness and kurtosis measures? Defaults to TRUE.
#' @return A bimodality coefficient (numeric)
#' @examples
#' data <- c(rnorm(15, 0, 1), rnorm(21, 5, 1))
#' hist(data)
#' bimodality_coefficient(data, TRUE)
bimodality_coefficient <- function(x, finite = TRUE)
{
  if (finite == TRUE) {
    G = UMtools::skewness(x, TRUE)
    K = UMtools::kurtosis(x, TRUE)
    n = length(x)
    B = ((G^2) + 1)/(K + ((3 * ((n - 1)^2))/((n - 2) * (n - 3))))
  }
  else {
    G = UMtools::skewness(x, FALSE)
    K = UMtools::kurtosis(x, FALSE)
    B = ((G^2) + 1)/K
  }
  return(B)
}
