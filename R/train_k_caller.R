#' K-caller performance for a given set of eps and minPts
#' @export
#' @import parallel
#' @import dbscan
#' @description It allows the quantification of the K-caller's performance for a given dataset, training set of probes
#' and parameters eps and minPts.
#' @param M Methylated fluorescence mean intensity matrix (CpGs as rows, samples as columns)
#' @param U Unmethylated fluorescence mean intensity matrix (CpGs as rows, samples as columns)
#' @param training_set List. As an example, see data(training_set). More details in help(training_set)
#' @param minPts dbscan parameter. Minimum numbers of points in the eps region to consider a core point; help(dbscan, dbscan) for more details.
#' @param eps dbscan parameter. Size of the neighbourhood; help(dbscan, dbscan) for more details.
#' @param nThread Number of CPU cores to employ
#' @return List. Containing confusion matrix and macro-Precission, macro-Recall and Macro F1-score for
#' 1-3 clusters and 1-4 clusters. Given the rarity of K=4 clusters, it is not fair to give the same weight
#' in evaluation of the multi-class classification performance. To decide what combination of parameters to employ,
#' we recommend to use macro F1-score for 1-3 clusters.
#' @examples
#' data(training_set)
#' train_k_caller(M, U, training_set, 3, 0.07)
#'
train_k_caller <- function(M, U, training_set, minPts, eps, nThread = NULL)
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

  k1_res = par_EW_Kcalling(M[training_set$k_1,], U[training_set$k_1,], minPts, eps, nThread)
  k2_res = par_EW_Kcalling(M[training_set$k_2,], U[training_set$k_2,], minPts, eps, nThread)
  k3_res = par_EW_Kcalling(M[training_set$k_3,], U[training_set$k_3,], minPts, eps, nThread)
  k4_res = par_EW_Kcalling(M[training_set$k_4,], U[training_set$k_4,], minPts, eps, nThread)
  pred = c(k1_res, k2_res, k3_res, k4_res)
  real = c(rep(1, length(training_set$k_1)), rep(2, length(training_set$k_2)),
           rep(3, length(training_set$k_3)), rep(4, length(training_set$k_4)))

  cm = table(real, pred)
  perfor1_3 = evaluate_multiclass(cm[1:3, 1:3]) # Please note that this excludes beyond K > 3 predictions.
  perfor1_4 = evaluate_multiclass(cm[1:4, 1:4]) # Please note that this excludes beyond K > 4 predictions.
  return(list(confusion_matrix = cm, performance_one_to_three_clusters = perfor1_3, performance_one_to_four_clusters = perfor1_4))
}
