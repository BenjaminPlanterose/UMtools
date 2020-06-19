#' Transforms G/R fluorescence signals into U/M
#'
#' @param Red A red fluorescence intensity matrix (probes as rows, samples as columns)
#' @param Grn A green fluorescence intensity matrix (probes as rows, samples as columns)
#' @param rgSet An rgSet object imported by minfi (see minfi::read.metharray.exp for details)
#' @return A list including matrices \code{M} and \code{U} (CpGs as rows, samples as columns)
#' @examples
#' GR_to_UM(Red, Grn, rgSet)
GR_to_UM <- function(Red, Grn, rgSet)
{
  locusNames <- getManifestInfo(rgSet, "locusNames")
  TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
  TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")
  TypeII <- getProbeInfo(rgSet, type = "II")
  M_and_U <- minfi:::.preprocessRaw(Red = Red, Green = Grn, locusNames = locusNames,
                                    TypeI.Red = TypeI.Red, TypeI.Green = TypeI.Green, TypeII = TypeII)
  list(M  = M_and_U[["M"]],
       U  = M_and_U[["U"]])
}
