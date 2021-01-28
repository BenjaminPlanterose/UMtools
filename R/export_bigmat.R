#' Export big matrices
#' @export
#' @import data.table
#' @description Exports large matrices with the help of data.table::fwrite
#' @param bigmat A matrix (probes as rows, samples as columns)
#' @param filename Name of the output file. The date will be appended automatically.
#' @param nThread Number of CPU cores to employ
#' @return NULL
#' @examples
#' export_bigmat(M, "M.txt")
export_bigmat <- function(bigmat, filename, nThread = 1)
{
  fwrite(data.table(bigmat, keep.rownames = T), paste(Sys.Date(), filename, sep = '_'), quote = F,
         row.names = F, col.names = T, sep = '\t', nThread = nThread)
}
