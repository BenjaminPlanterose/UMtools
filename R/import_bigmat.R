#' Import big matrices
#' @export
#' @import data.table
#' @description Imports large matrices with the help of data.table::fread
#' @param filename Name of the file to import
#' @param nThread Number of CPU cores to employ
#' @return Matrix
#' @examples
#' M = import_bigmat("2021-01-27_M.txt")
import_bigmat <- function(filename, nThread = 1)
{
  bigmat = fread(filename, nThread = nThread)
  cpgs <- as.vector(bigmat$rn)
  bigmat <- bigmat[, !"rn"]
  bigmat <- as.matrix(bigmat)
  rownames(bigmat) <- cpgs
  return(bigmat)
}
