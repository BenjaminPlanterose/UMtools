#' Computes bimodality per CpG across samples in a coefficient of variation matrix
#' @description a
#' @details a
#' @param CV A coefficient of variation matrix (probes as rows, samples as columns)
#' @return A vector of coefficient of bimodality per CpG, \code{BC}
#' @examples
#' compute_BC_CV(CV)
par_EW_Kcalling <- function(M, U, minPts = 25, reach = seq(0.328, 0.352, 0.004), R = 2)
{
  np <- detectCores(logical = FALSE)
  cl <- makeCluster(np)
  clusterExport(cl, c("k_per_CpG", "dbscan", "R", "reach", "minPts"), envir=environment())

  print(dim(M))
  print(dim(U))

  r <- parSapply(cl = cl, X = 1:nrow(M), FUN = function(X) k_per_CpG(data.frame(x = M[X,], y = U[X,]),
                                                                     minPts, reach, R)) # process per row
  names(r) <- rownames(M)

  stopCluster(cl)
  return(r)
}
