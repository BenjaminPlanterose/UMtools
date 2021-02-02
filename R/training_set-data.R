#' K-caller (Training Set; EUR, Illumina Infinium HumanMethylation450 Beadchip microarray)
#' @docType data
#' @usage data(training_set)
#' @format An object of class list
#' @keywords datasets
#' @description  Includes a list 4 vectors of 450K cg id's:
#' \itemize{
#'   \item K = 1
#'   \item K = 2
#'   \item K = 3
#'   \item K = 4
#' }
#' These correspond to probes that form from 1 to 4 clusters in the U/M plane on a dataset of
#' 426 MZ twin pairs of the E-risk cohort (British, EUR ancestry). These sets were compiled
#' by manually examining U/M plots of random probes. However, the number of clusters observed
#' may be different in other datasets. Key parameters that will determine if this is true
#'  will be  sample size or ancestry.
#' If one manually verifies that these lists apply to their own dataset, they can be used for
#' calibrating parameters \emph{eps} and \emph{minPts} from the K-caller.
#' @references Revisiting Genetic artefacts on DNA methylation microarrays. Genome Research
#' @examples
#' data(training_set)
"training_set"
