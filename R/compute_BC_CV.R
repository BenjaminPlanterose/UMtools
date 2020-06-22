#' Computes bimodality per CpG across samples in a coefficient of variation matrix
#' @description a
#' @details a
#' @param CV A coefficient of variation matrix (probes as rows, samples as columns)
#' @return A vector of coefficient of bimodality per CpG, \code{BC}
#' @examples
#' compute_BC_CV(CV)
compute_BC_CV <- function(CV)
{
  BC = sapply(1:nrow(CV), function(x) bimodality_coefficient(CV[x, ]))
  names(BC) = rownames(CV)
  return(BC)
}
