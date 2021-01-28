#' K-caller
#' @export
#' @import dbscan
#' @import RColorBrewer
#' @description Employing preprocessing in combination with the algorithm dbscan, it allows the counting of
#' clusters in the U/M plane. If plot = TRUE, a U/M plot is produced for which each cluster in coloured. Outliers,
#' samples belonging to no cluster, will be coloured in black (as defined by the dbscan algorithm).
#' @details Prior tuning of parameters {eps, minPts} is absolutely necessary and will depend on the sample size of the dataset employed.
#' Parameters minPts = 12, eps = 0.05 were optimal in the E-risk cohort (n = 852). Smaller datasets will probably require
#' reducing minPts and increasing eps, while larger datasets will probably require increasing minPts and reducing eps.
#' For finding optimal parameters, are training set is available at UMtools; see help(training_set) for more details.
#' @param CpG Targeted CpG
#' @param M Methylated fluorescence mean intensity matrix (CpGs as rows, samples as columns)
#' @param U Unmethylated fluorescence mean intensity matrix (CpGs as rows, samples as columns)
#' @param plot Logical. Whether to visualize the K-callers output in the U/M plane or not.
#' @param minPts dbscan parameter. Minimum numbers of points in the eps region to consider a core point; help(dbscan, dbscan) for more details.
#' @param eps dbscan parameter. Size of the neighbourhood; help(dbscan, dbscan) for more details.
#' @return Numeric, number of predicted clusters in the U/M plane.
#' @examples
#' Kcall_CpG("cg03398919", M, U)
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
    plot(m, u, col = as.character(col), pch = 19, main = CpG, xlab = "M", ylab = "U")
  }
  nclust = length(table(res[res != 0]))
  return(nclust)
}

