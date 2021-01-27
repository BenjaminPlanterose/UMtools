#' Computes coefficient of variation of the log of the total intensity signal
#' @description a
#' @details a
#' @param M_SD Methylated standard deviation intensity matrix (CpGs as rows, samples as columns)
#' @param U_SD Unmethylated fluorescence standard deviation intensity matrix (CpGs as rows, samples as columns)
#' @param M Methylated fluorescence mean intensity matrix (CpGs as rows, samples as columns)
#' @param U Unmethylated fluorescence mean intensity matrix (CpGs as rows, samples as columns)
#' @return A matrix \code{CV} (CpGs as rows, samples as columns)
#' @examples
#' rgSet = read.metharray.exp(getwd(), extended = TRUE)
#' Grn = assay(rgSet, "Green")       # Green mean across beads
#' Red = assay(rgSet, "Red")         # Red mean across beads
#' M_U = GR_to_UM(Red, Grn, rgSet)
#' GrnSD = assay(rgSet, "GreenSD")   # Green SD across beads
#' RedSD = assay(rgSet, "RedSD")     # Red SD across beads
#' GR_to_UM(RedSD, GrnSD, rgSet)
#' compute_CV(M_U_sd$M, M_U_sd$U, M_U$M, M_U_sd$U)
compute_CV = function(M_SD, U_SD, M, U, alpha = 100)
{
  R = (M_SD + U_SD + alpha)/(U + M + alpha)
  1/(log(U + M + alpha)/R - R/2)
}
