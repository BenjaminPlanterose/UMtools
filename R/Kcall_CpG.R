#' Transforms number of beads matrix from G/R to U/M
#' @description a
#' @details a
#' @param nBeads A matrix containing the number of beads per probe and per sample (probes as rows, samples as columns)
#' @param rgSet An rgSet object imported by minfi (see minfi::read.metharray.exp for details)
#' @return A matrix containing number of beads per CpG and per sample, \code{nBeads_cg}, (CpGs as rows, samples as columns)
#' @examples
#' beads_GR_to_UM(nBeads, rgSet)
Kcall_CpG <- function(CpG, M, U, minPts = 25, plot = T, reach = seq(0.328, 0.352, 0.004))
{
  m = M[CpG,]; u = U[CpG,]
  df = data.frame(x = 1 - atan((u+100)/(m+100)), y = log2(u + m + 100))
  df = as.data.frame(scale(df, center = T, scale = T))
  if(plot)
  {
    db <- dbscan::dbscan(x = df, eps = mean(reach), minPts)
    plot(df$x, df$y, col = as.factor(db$cluster), pch = 19, main = CpG)
  }

  res = lapply(X = reach, FUN = function(x) dbscan::dbscan(df, eps = x, minPts)$cluster)
  res = lapply(res, table)
  res = lapply(res, function(x) x[names(x) != "0"])
  res = sapply(res, length)
  nclust = min(as.numeric(res))
  return(nclust)
}
