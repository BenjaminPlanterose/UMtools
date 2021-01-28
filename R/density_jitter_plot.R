#' Density-jitter plot
#' @export
#' @import modes
#' @description Overlaps a density plot with a jitter plot. A random y-component is added to cause the spread of
#' points.
#' @param mat Matrix of a continous varible such as beta-value, M-value, CVlogT, etc (CpGs as rows and samples as columns)
#' @param CpG Targeted CpG
#' @param sex vector of sex. First level is expected to be Female (coloured pink). Check levels(sex) to see if that is the case.
#' @param alpha Transparency
#' @param lambda Regulates the random spread of the points
#' @param main Title name

#' @return Graphics
#' @examples
#' density_jitter_plot(beta_value, "cg00050873", pheno$sex)
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
