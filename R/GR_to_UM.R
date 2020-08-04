#' GR_to_UM
#' @description Transforms Green/Red fluorescence signals into Unmethylated/Methylated channels
#' @details
#'Illumina DNA methylation microarrays detection rely on three types of probes:
#' \itemize{
#'   \item type-I Red: 2 probes/addresses per CpG, informative on the Red channel exclusively.
#'   \item type-I Green: 2 probes/addresses per CpG, informative on the Green channel exclusively.
#'   \item type-II: 1 probe/address per CpG, informative in both Green and Red channels.
#' }
#'
#' When *what* is equal to "SD" or "Mean", Methylated and unmethylated intensities are computed as followed:
#'
#' |                   | **Methylated**   | **Unmethylated** |
#' |:----------------- |:----------------:| ----------------:|
#' | **Type-II**       | Green (addressA) | Red (addressA)   |
#' | **Type-I Red**    | Red (addressB)   | Red (addressA)   |
#' | **Type-I Green**  | Green (addressB) | Green (addressA) |
#'
#'Also, due to the beadchip technology, each probe has been assigned to a random number of beads.
#'If *what* == "NBeads", it converts the number of beads per probe and per sample to
#'the number of beads per CpG and per sample. This requires a criteria for type-I probes for
#'which 2 probes target the same CpG. In this case,the minimum number between the two probes
#'targetting a type-I probe is chosen.
#'
#' @md
#' @param Red Mean or SD Red fluorescence intensity matrix (probes as rows, samples as columns)
#' @param Grn Mean or SD Green fluorescence intensity matrix (probes as rows, samples as columns)
#' @param nBeads A matrix containing the number of beads (probes as rows, samples as columns)
#' @param rgSet An rgSet object imported by minfi (see minfi::read.metharray.exp for details)
#' @param what Either 'Mean', 'SD' or 'NBeads'
#' @return A list including matrices *M* and *U* (CpGs as rows, samples as columns) when *what* is "Mean"
#' or "SD" and or a matrix *nBeads*, containing number of beads (CpGs as rows, samples as columns) when
#' *what* is equal to "NBeads".
#' @examples
#' # Read IDAT files. Without the extended option, SD cannot be acquired
#' rgSet <- read.metharray.exp(getwd(), extended = TRUE)
#'
#' # Green/Red mean to Unmethylated/Methylated mean intensities
#' Grn <- assay(rgSet, "Green")
#' Red <- assay(rgSet, "Red")
#' UM <- GR_to_UM(Red = Red, Grn = Grn, rgSet = rgSet, what = "Mean")
#'
#' # Green/Red SD to Unmethylated/Methylated intensity SDs
#' GrnSD <- assay(rgSet, "GreenSD")
#' RedSD <-assay(rgSet, "RedSD")
#' UM_SD <- GR_to_UM(Red = RedSD, Grn = GrnSD, rgSet = rgSet, what = "SD")
#'
#' # nBeads in probes to nBeads per CpG
#' nBeads <- assay(rgSet, "NBeads")
#' nBeads <- GR_to_UM(nBeads = nBeads, rgSet = rgSet, what = "NBeads")
GR_to_UM <- function(Red = NULL, Grn = NULL, rgSet, what = NULL, nBeads = NULL)
{
  if(is.null(what)) stop("'what' needs to be equal to 'SD', 'Mean' or 'NBeads'. Please see documentation")
  locusNames <- minfi::getManifestInfo(rgSet, "locusNames")
  TypeI.Red <- minfi::getProbeInfo(rgSet, type = "I-Red")
  TypeI.Green <- minfi::getProbeInfo(rgSet, type = "I-Green")
  TypeII <- minfi::getProbeInfo(rgSet, type = "II")

  if(what %in% c("Mean", "SD"))
  {
    M_and_U <- minfi:::.preprocessRaw(Red = Red, Green = Grn, locusNames = locusNames,
                                      TypeI.Red = TypeI.Red, TypeI.Green = TypeI.Green, TypeII = TypeII)
    list(M = M_and_U[["M"]],
         U = M_and_U[["U"]])
  }
  else if(what == "NBeads")
  {
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
}
