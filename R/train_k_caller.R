#' Transforms number of beads matrix from G/R to U/M
#' @description a
#' @details a
#' @param nBeads A matrix containing the number of beads per probe and per sample (probes as rows, samples as columns)
#' @param rgSet An rgSet object imported by minfi (see minfi::read.metharray.exp for details)
#' @return A matrix containing number of beads per CpG and per sample, \code{nBeads_cg}, (CpGs as rows, samples as columns)
#' @examples
#' beads_GR_to_UM(nBeads, rgSet)
#'
train_k_caller <- function(M, U, training_set, minPts, eps)
{
  evaluate_multiclass <- function(cm) # Obtained from the following post: https://blog.revolutionanalytics.com/2016/03/com_class_eval_metrics_r.html
  {
    n = sum(cm) # number of instances
    nc = nrow(cm) # number of classes
    diag = diag(cm) # number of correctly classified instances per class
    rowsums = apply(cm, 1, sum) # number of instances per class
    colsums = apply(cm, 2, sum) # number of predictions per class
    p = rowsums / n # distribution of instances over the actual classes
    q = colsums / n # distribution of instances over the predicted
    accuracy = sum(diag) / n
    precision = diag / colsums
    recall = diag / rowsums
    f1 = 2 * precision * recall / (precision + recall)
    macroPrecision = mean(precision)
    macroRecall = mean(recall)
    macroF1 = mean(f1)

    return(c(macroPrecision = macroPrecision, macroRecall = macroRecall, macroF1 = macroF1))
  }

  k1_res = par_EW_Kcalling(M[training_set$k_1,], U[training_set$k_1,], minPts, eps)
  k2_res = par_EW_Kcalling(M[training_set$k_2,], U[training_set$k_2,], minPts, eps)
  k3_res = par_EW_Kcalling(M[training_set$k_3,], U[training_set$k_3,], minPts, eps)
  k4_res = par_EW_Kcalling(M[training_set$k_4,], U[training_set$k_4,], minPts, eps)
  pred = c(k1_res, k2_res, k3_res, k4_res)
  real = c(rep(1, length(k_1)), rep(2, length(k_2)),
           rep(3, length(k_3)), rep(4, length(k_4)))

  cm = table(real, pred)
  perfor1_3 = evaluate_multiclass(cm[1:3, 1:3])
  perfor1_4 = evaluate_multiclass(cm[1:4, 1:4])
  return(list(confusion_matrix = cm, performance_one_to_three_clusters = perfor1_3, performance_one_to_four_clusters = perfor1_4))
}
