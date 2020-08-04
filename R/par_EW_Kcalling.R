#' Computes bimodality per CpG across samples in a coefficient of variation matrix
#' @description a
#' @details a
#' @param CV A coefficient of variation matrix (probes as rows, samples as columns)
#' @return A vector of coefficient of bimodality per CpG, \code{BC}
#' @examples
#' compute_BC_CV(CV)
par_EW_Kcalling <- function(M, U, minPts = 12, eps = 0.035)
{
  k_per_CpG <- function(df, minPts = 12, eps = 0.035)
  {
    m = df$x; u = df$y
    df = data.frame(x = (m + 0)/(u + m + 100), y = log2(u + m + 100))
    df$y = (df$y - min(df$y))/(max(df$y)) # Switch
    res = lapply(X = eps, FUN = function(x) dbscan::dbscan(df, eps = x, minPts)$cluster)
    res = lapply(res, table)
    res = lapply(res, function(x) x[names(x) != "0"])
    res = sapply(res, length)
    nclust = min(as.numeric(res))
    return(nclust)
  }

  np <- detectCores(logical = FALSE)
  cl <- makeCluster(np)
  clusterExport(cl, c("k_per_CpG", "dbscan", "eps", "minPts"), envir=environment())
  print(dim(M)); print(dim(U))
  r <- parSapply(cl = cl, X = 1:nrow(M), FUN = function(X) k_per_CpG(data.frame(x = M[X,], y = U[X,]),
                                                                     minPts, eps)) # process per row
  names(r) <- rownames(M)
  stopCluster(cl)
  return(r)
}
