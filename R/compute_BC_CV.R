#' Computes bimodality per CpG across samples in a coefficient of variation matrix
#' @description a
#' @details a
#' @param CV A coefficient of variation matrix (probes as rows, samples as columns)
#' @return A vector of coefficient of bimodality per CpG, \code{BC}
#' @examples
#' compute_BC_CV(CV)
compute_BC_CV <- function(CV, ncores = NULL)
{
  if(is.null(ncores))
  {
    np <- detectCores(logical = FALSE)
  }

  cl <- makeCluster(np)
  clusterExport(cl, c("bimodality_coefficient"), envir=environment())
  r <- parSapply(cl = cl, X = 1:nrow(CV), FUN = function(x) bimodality_coefficient(CV[x, ]))
  stopCluster(cl)

  names(r) <- rownames(CV)
  return(r)
}

