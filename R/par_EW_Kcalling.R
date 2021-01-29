#' Epigenome-wide K-calling
#' @export
#' @import dbscan
#' @import parallel
#' @description Same as Kcall_CpG but it performs K-calling on all the CpGs (rows) of a given matrix.
#' Also, it allows parallel processing. If nThread is not specified, it will automatically detect
#' how many clusters are available in the machine and will use N-1 cores.
#' @param M Methylated fluorescence mean intensity matrix (CpGs as rows, samples as columns)
#' @param U Unmethylated fluorescence mean intensity matrix (CpGs as rows, samples as columns)
#' @param minPts dbscan parameter. Minimum numbers of points in the eps region to consider a core point; help(dbscan, dbscan) for more details.
#' @param eps dbscan parameter. Size of the neighbourhood; help(dbscan, dbscan) for more details.
#' @param nThread Number of CPU cores to employ
#' @return Vector. Predicted number of clusters for each CpG.
#' @examples
#' annotation <- getAnnotation(rgSet)
#' chrY = rownames(annotation)[annotation$chr == "chrY"]
#' K_vec = par_EW_Kcalling(M[chrY,], U[chrY,])
#'
par_EW_Kcalling <- function(M, U, minPts = 12, eps = 0.035, nThread = NULL)
{

  k_per_CpG <- function(df, minPts = 12, eps = 0.035)
  {
    m = df$x; u = df$y
    df = data.frame(x = m/(u + m + 100), y = log2(u + m + 100))
    df$y = (df$y - min(df$y))/(max(df$y)) # Switch
    res = dbscan::dbscan(df, eps, minPts)$cluster
    nclust = length(table(res[res != 0]))
    return(nclust)
  }

  if(is.null(nThread))
  {
    nThread <- parallel::detectCores(logical = FALSE) - 1
  }

  print(dim(M)); print(dim(U)); message(paste("Using", nThread, "cores"))
  r <- parallel::mclapply(1:nrow(M), function(X) k_per_CpG(data.frame(x = M[X,], y = U[X,]), minPts, eps),
                          mc.cores = nThread)
  names(r) <- rownames(M)
  return(unlist(r))
}
