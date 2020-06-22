#' Transforms number of beads matrix from G/R to U/M
#' @description a
#' @details a
#' @param nBeads A matrix containing the number of beads per probe and per sample (probes as rows, samples as columns)
#' @param rgSet An rgSet object imported by minfi (see minfi::read.metharray.exp for details)
#' @return A matrix containing number of beads per CpG and per sample, \code{nBeads_cg}, (CpGs as rows, samples as columns)
#' @examples
#' beads_GR_to_UM(nBeads, rgSet)
export_bigmat <- function(bigmat, filename, nThread = 1)
{
  fwrite(data.table(bigmat, keep.rownames = T), paste(Sys.Date(), filename, sep = '_'), quote = F,
         row.names = F, col.names = T, sep = '\t', nThread = nThread)
}
