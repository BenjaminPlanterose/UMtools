#' Transforms G/R fluorescence signals into U/M
#' @description a
#' @details a
#' @param Red A red fluorescence intensity matrix (probes as rows, samples as columns)
#' @param Grn A green fluorescence intensity matrix (probes as rows, samples as columns)
#' @param rgSet An rgSet object imported by minfi (see minfi::read.metharray.exp for details)
#' @return A list including matrices \code{M} and \code{U} (CpGs as rows, samples as columns)
#' @examples
#' GR_to_UM(Red, Grn, rgSet)
UM_plot <- function(M, U, CpG, sex = NULL, xlim = NULL, ylim = NULL, alpha = 0.5)
{
  df = data.frame(x = M[CpG,], y = U[CpG,])

  if(is.null(xlim))
  {
    thr = max(c(df$x, df$y))
    nDig = floor(log10(thr)) + 1
    xlim = c(0, 10^(nDig-2)*ceiling(thr/10^(nDig-2)))
  }

  if(is.null(ylim))
  {
    thr = max(c(df$x, df$y))
    nDig = floor(log10(thr)) + 1
    ylim = c(0, 10^(nDig-2)*ceiling(thr/10^(nDig-2)))
  }

  if(is.null(sex))
  {
    plot(df$x, df$y, col = alpha("black", alpha), pch = 19, main = CpG, xlab = 'M', ylab = 'U',
         xlim = xlim, ylim = ylim, cex.main = 0.8)
  }
  else
  {
    col = as.factor(sex)
    levels(col) = c("deeppink2", "dodgerblue1")
    col = as.character(col)
    plot(df$x, df$y, col = alpha(col, alpha), pch = 19, main = CpG, xlab = 'M', ylab = 'U',
         xlim = xlim, ylim = ylim, cex.main = 0.8)
  }

}
