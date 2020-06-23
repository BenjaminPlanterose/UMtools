#' Computes bimodality per CpG across samples in a coefficient of variation matrix
#' @description a
#' @details a
#' @param CV A coefficient of variation matrix (probes as rows, samples as columns)
#' @return A vector of coefficient of bimodality per CpG, \code{BC}
#' @examples
#' compute_BC_CV(CV)
k_per_CpG <- function(df, minPts = 25, reach = seq(0.328, 0.352, 0.004), R = 2)
{
  x = df$x
  y = df$y
  alpha1 = IQR(x)/R
  alpha2 = IQR(y)/R

  df = data.frame(x = log2(x + alpha1), y = log2(y + alpha2))
  df = as.data.frame(scale(df, center = T, scale = T))

  res = lapply(X = reach, FUN = function(x) dbscan::dbscan(df, eps = x, minPts)$cluster)
  res = lapply(res, table)
  res = lapply(res, function(x) x[names(x) != "0"])
  res = sapply(res, length)
  nclust = min(as.numeric(res))

  return(nclust)
}
