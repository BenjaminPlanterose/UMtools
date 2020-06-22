#' Transforms G/R fluorescence signals into U/M
#' @description a
#' @details a
#' @param Red A red fluorescence intensity matrix (probes as rows, samples as columns)
#' @param Grn A green fluorescence intensity matrix (probes as rows, samples as columns)
#' @param rgSet An rgSet object imported by minfi (see minfi::read.metharray.exp for details)
#' @return A list including matrices \code{M} and \code{U} (CpGs as rows, samples as columns)
#' @examples
#' GR_to_UM(Red, Grn, rgSet)
import_bigmat <- function(filename, nThread = 1)
{
  bigmat = fread(filename, nThread = nThread)
  cpgs <- as.vector(bigmat$rn)
  bigmat <- bigmat[, !"rn"]
  bigmat <- as.matrix(bigmat)
  rownames(bigmat) <- cpgs
  return(bigmat)
}
