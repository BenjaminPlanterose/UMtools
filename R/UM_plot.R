#' U/M plot
#' @export
#' @import scales
#' @description Scatter plot of fluorescence intensities: unmethylated (U) fluorescence against methylated (M) fluorescence.
#' @details More information on what it is referred as U/M fluorescence is available at help(GR_to_UM).
#' @param M Methylated fluorescence mean intensity matrix (CpGs as rows, samples as columns)
#' @param U Unmethylated fluorescence mean intensity matrix (CpGs as rows, samples as columns)
#' @param CpG A CpG identifier (e.g. "cg15771735")
#' @param sex vector of sex. First level is expected to be Female (coloured pink). Check levels(sex) to see if that is the case.
#' @param xlim x-axis limit
#' @param ylim y-axis limit
#' @param alpha transparency
#' @return Graphics
#' @examples
#' UM_plot(M = M, U = U, CpG = "cg00050873", sex = pheno$sex)
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
         xlim = xlim, ylim = ylim)
  }
  else
  {
    col = as.factor(sex)
    levels(col) = c("deeppink2", "dodgerblue1")
    col = as.character(col)
    plot(df$x, df$y, col = alpha(col, alpha), pch = 19, main = CpG, xlab = 'M', ylab = 'U',
         xlim = xlim, ylim = ylim)
  }

}
