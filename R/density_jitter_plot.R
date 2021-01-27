#' Plots density jitter plot
#' @description a
#' @details a
#' @param Red A red fluorescence intensity matrix (probes as rows, samples as columns)
#' @param Grn A green fluorescence intensity matrix (probes as rows, samples as columns)
#' @param rgSet An rgSet object imported by minfi (see minfi::read.metharray.exp for details)
#' @return A list including matrices \code{M} and \code{U} (CpGs as rows, samples as columns)
#' @examples
#' GR_to_UM(Red, Grn, rgSet)
density_jitter_plot <- function(mat, CpG, sex = NULL, alpha = 0.5, lambda = 0.05, main = NULL)
{
  dens = density(mat[CpG,])
  mu = median(dens$y)
  sd = diff(quantile(dens$y, c(0.5-lambda, 0.5 + lambda)))

  if(is.null(sex))
  {
    plot(dens, main = CpG, xlab = deparse(substitute(mat)))
    points(x = mat[CpG,], y = rnorm(ncol(mat), mu, sd), pch = 19)
  }

  else if(!is.null(sex))
  {
    col = as.factor(sex)
    levels(col) = c("deeppink2", "dodgerblue1")
    col = as.character(col)
    plot(dens, xlab = deparse(substitute(mat)), main = CpG)
    points(x = mat[CpG,], y = rnorm(ncol(mat), mu, sd), col = alpha(col, alpha), pch = 19)
  }

}
