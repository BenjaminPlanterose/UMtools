#' BC(CV)
#' @import parallel
#' @export
#' @description Computes bimodality coefficient of the CVlogT per CpG across samples
#' @details It computes the following per row (per CpG):
#' \deqn{BC = \frac{\gamma^{2} +1}{\kappa + \frac{3\cdot(n-1)^{2}}{(n-2)\cdot(n-3)}}}
#' Where \eqn{\gamma} is the sample skewness and where \eqn{\kappa} is the sample excess kurtosis.
#'
#' For more details, check the reference paper.
#' @references Revisiting Genetic artefacts on DNA methylation microarrays. Genome Research
#' @param CV A coefficient of variation matrix (probes as rows, samples as columns)
#' @param parallel Whether to perform this task employing parallel processing.
#' @param nThread Number of CPU cores to employ on this task
#' @return A vector of coefficients of bimodality per CpG
#' @examples
#' rgSet = read.metharray.exp(getwd(), extended = TRUE)
#' Grn = assay(rgSet, "Green")       # Green mean across beads
#' Red = assay(rgSet, "Red")         # Red mean across beads
#' M_U = GR_to_UM(Red, Grn, rgSet, "Mean")
#' GrnSD = assay(rgSet, "GreenSD")   # Green SD across beads
#' RedSD = assay(rgSet, "RedSD")     # Red SD across beads
#' M_U_sd = GR_to_UM(RedSD, GrnSD, rgSet, "SD")
#' CV = compute_CV(M_U_sd$M, M_U_sd$U, M_U$M, M_U_sd$U)
#' compute_BC_CV(CV)
compute_BC_CV <- function(CV, parallel = F, nThread = NULL)
{
  if(parallel == T)
  {
    if(is.null(nThread))
    {
      nThread <- detectCores(logical = FALSE)
    }
    print(dim(CV)); message(paste("Using", nThread, "cores"))
    r <- parallel::mclapply(1:nrow(CV), function(x) UMtools::bimodality_coefficient(CV[x, ]), mc.cores = nThread)
    r <- unlist(r)
  }
  else
  {
    r = sapply(1:nrow(CV), function(x) UMtools::bimodality_coefficient(CV[x, ]))
  }

  names(r) <- rownames(CV)
  return(r)
}

