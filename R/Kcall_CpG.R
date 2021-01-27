#' Transforms number of beads matrix from G/R to U/M
#' @description a
#' @details a
#' @param nBeads A matrix containing the number of beads per probe and per sample (probes as rows, samples as columns)
#' @param rgSet An rgSet object imported by minfi (see minfi::read.metharray.exp for details)
#' @return A matrix containing number of beads per CpG and per sample, \code{nBeads_cg}, (CpGs as rows, samples as columns)
#' @examples
#' beads_GR_to_UM(nBeads, rgSet)
Kcall_CpG <- function(CpG, M, U, plot = T, minPts = 12, eps = 0.035)
{
  m = M[CpG,]; u = U[CpG,]
  df = data.frame(x = m/(u + m + 100), y = log2(u + m + 100))
  df$y = (df$y - min(df$y))/(max(df$y)) # Switch
  res = dbscan::dbscan(df, eps, minPts)$cluster
  if(plot)
  {
    K = length(table(res[res != 0]))
    if(K > 8)
    {
      stop("cannot plot; more clusters than colours")
    }
    pal = brewer.pal(8, "Dark2")[1:K]
    col = as.factor(res)
    if(sum(col == "0") > 0)
    {
      levels(col) = c("black", pal)
    }
    else
    {
      levels(col) = pal
    }
    plot(m, u, col = as.character(col), pch = 19, main = CpG)
  }
  nclust = length(table(res[res != 0]))
  return(nclust)
}
