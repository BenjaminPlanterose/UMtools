#' Transforms the matrix of number of beads from Green/Red to Unmethylated/Methylated channels
#' @description Transforms
#' @details
#'Illumina DNA methylation microarrays detection rely on three types of probes:
#' \itemize{
#'   \item type-I Red: 2 probes/addresses per CpG, informative on the Red channel exclusively.
#'   \item type-I Green: 2 probes/addresses per CpG, informative on the Green channel exclusively.
#'   \item type-II: 1 probe/address per CpG, informative in both channels.
#' }
#'Due to the beadchip technology, each probe has been assigned to a random number of beads.
#'beads_GR_to_UM converts the number of beads per probe and per sample to
#'the number of beads per CpG and per sample.
#'This requires a criteria for type-I probes for which 2 probes target the same CpG. In this case,
#'the minimum number between the two probes targetting a type-I probe is chosen.
#'
#'
#' @param nBeads A matrix containing the number of beads per probe and per sample (probes as rows, samples as columns)
#' @param rgSet An rgSet object imported by minfi (see minfi::read.metharray.exp for details)
#' @return A matrix containing number of beads per CpG and per sample, \code{nBeads_cg}, (CpGs as rows, samples as columns)
#' @examples
#' beads_GR_to_UM(nBeads, rgSet)
beads_GR_to_UM <- function(nBeads, rgSet)
{
  locusNames <- getManifestInfo(rgSet, "locusNames")
  TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
  TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")
  TypeII <- getProbeInfo(rgSet, type = "II")

  nBeads_cg = matrix(NA, nrow = length(locusNames), ncol = ncol(rgSet))
  rownames(nBeads_cg) = locusNames
  colnames(nBeads_cg) = colnames(rgSet)
  nBeads_cg[TypeII$Name,] = nBeads[TypeII$AddressA,]
  comp = nBeads[TypeI.Red$AddressA,] < nBeads[TypeI.Red$AddressB,]
  nBeads_cg[TypeI.Red$Name,] = comp*nBeads[TypeI.Red$AddressA,] + (!comp)*nBeads[TypeI.Red$AddressB,]
  comp = nBeads[TypeI.Green$AddressA,] < nBeads[TypeI.Green$AddressB,]
  nBeads_cg[TypeI.Green$Name,] = comp*nBeads[TypeI.Green$AddressA,] + (!comp)*nBeads[TypeI.Green$AddressB,]
  return(nBeads_cg)
}

