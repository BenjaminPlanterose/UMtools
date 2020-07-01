#' Computes BC(CV)
#' @description Computes bimodality coefficient of the coefficient of variation per CpG across samples
#' @details a
#' @param CV A coefficient of variation matrix (probes as rows, samples as columns)
#' @param nCores Number of cores to employ on this task
#' @return A vector of coefficient of bimodality per CpG, \code{BC}
#' @examples
#' rgSet = read.metharray.exp(getwd(), extended = TRUE)
#' Grn = assay(rgSet, "Green")       # Green mean across beads
#' Red = assay(rgSet, "Red")         # Red mean across beads
#' M_U = GR_to_UM(Red, Grn, rgSet)
#' GrnSD = assay(rgSet, "GreenSD")   # Green SD across beads
#' RedSD = assay(rgSet, "RedSD")     # Red SD across beads
#' GR_to_UM(RedSD, GrnSD, rgSet)
#' CV = compute_CV(M_U_sd$M, M_U_sd$U, M_U$M, M_U_sd$U)
#' compute_BC_CV(CV)
compute_BC_CV <- function(CV, nCores = NULL)
{
  if(is.null(nCores))
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

