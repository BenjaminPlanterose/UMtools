#' Transforms Green/Red fluorescence signals into Unmethylated/Methylated channels
#' @description Transforms
#' @details
#'Illumina DNA methylation microarrays detection rely on three types of probes:
#' \itemize{
#'   \item type-I Red: 2 probes/addresses per CpG, informative on the Red channel exclusively.
#'   \item type-I Green: 2 probes/addresses per CpG, informative on the Green channel exclusively.
#'   \item type-II: 1 probe/address per CpG, informative in both Green and Red channels.
#' }
#'
#' Methylated and unmethylated intensities are computed as followed:
#'
#' |                   | **Methylated**   | **Unmethylated** |
#' |:----------------- |:----------------:| ----------------:|
#' | **Type-II**       | Green (addressA) | Red (addressA)   |
#' | **Type-I Red**    | Red (addressB)   | Red (addressA)   |
#' | **Type-I Green**  | Green (addressB) | Green (addressA) |
#'
#' @md
#' @param Red A red fluorescence intensity matrix (probes as rows, samples as columns)
#' @param Grn A green fluorescence intensity matrix (probes as rows, samples as columns)
#' @param rgSet An rgSet object imported by minfi (see minfi::read.metharray.exp for details)
#' @return A list including matrices \code{M} and \code{U} (CpGs as rows, samples as columns)
#' @examples
#' rgSet = read.metharray.exp(getwd(), extended = TRUE)
#' Grn = assay(rgSet, "Green")       # Green mean across beads
#' Red = assay(rgSet, "Red")         # Red mean across beads
#' GR_to_UM(Red, Grn, rgSet)
#' GrnSD = assay(rgSet, "GreenSD")   # Green SD across beads
#' RedSD = assay(rgSet, "RedSD")     # Red SD across beads
#' GR_to_UM(RedSD, GrnSD, rgSet)
GR_to_UM <- function(Red, Grn, rgSet)
{
  locusNames <- minfi::getManifestInfo(rgSet, "locusNames")
  TypeI.Red <- minfi::getProbeInfo(rgSet, type = "I-Red")
  TypeI.Green <- minfi::getProbeInfo(rgSet, type = "I-Green")
  TypeII <- minfi::getProbeInfo(rgSet, type = "II")
  M_and_U <- minfi:::.preprocessRaw(Red = Red, Green = Grn, locusNames = locusNames,
                                    TypeI.Red = TypeI.Red, TypeI.Green = TypeI.Green, TypeII = TypeII)
  list(M  = M_and_U[["M"]],
       U  = M_and_U[["U"]])
}
