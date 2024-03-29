#' Bivariate Gaussian Mixture Model
#' @export
#' @import EMCluster
#' @import scales
#' @import RColorBrewer
#' @description Fits a Gaussian Mixture Model (bGMM) with K components on the Unmethylated/Methylated plane
#' @details It employs the function init.EM from the R-package EMCluster that fits a Gaussian Mixture Model employing the expectation-maximization
#' algorithm on the U/M plane. Run help(init.EM) for more information on the EMCluster control functions.
#' @param M A matrix containing the methylated intensities (CpGs as rows, samples as columns)
#' @param U A matrix containing the unmethylated intensities (CpGs as rows, samples as columns)
#' @param CpG A CpG identifier (e.g. "cg15771735")
#' @param K Targeted number of clusters
#' @param ... Control functions for EMCluster
#' @return A vector of classes in the same order as the columns of U and M.
#' @examples
#' rgSet = read.metharray.exp(getwd(), extended = TRUE)
#' Grn = assay(rgSet, "Green")       # Green mean across beads
#' Red = assay(rgSet, "Red")         # Red mean across beads
#' M_U = GR_to_UM(Red, Grn, rgSet)
#' bGMM(M = M_U$M, U = M_U$U, CpG = "cg00814218", K = 3)
bGMM <- function(M, U, CpG, K, stable.solution = TRUE, min.n = NULL, min.n.iter = 2000, method = 'em.EM', EMC = .EMC, transform = TRUE)
{
  m = M[CpG, ]; m = m + rnorm(n = length(m), mean = 0, sd = 0.1)
  u = U[CpG, ]; u = u + rnorm(n = length(u), mean = 0, sd = 0.1)
  
  if (transform) 
  {
    df = data.frame(x = m/(u + m + 100), y = log2(u + m +  100))
    df$y = (df$y - min(df$y))/(max(df$y))
  }
  else
  {
    df = data.frame(x = m, y = u)
  }
  
  ret <- EMCluster::init.EM(x = df, nclass = K, method = "em.EM", 
                            min.n.iter = min.n.iter, lab = NULL, EMC = .EMC, stable.solution = stable.solution, 
                            min.n = min.n)
  class_i <- EMCluster::assign.class(df, ret, return.all = FALSE)$class
  if (K == 1 | K > 4) {
    if (K > 8) {
      stop("Too many K and not enough colours")
    }
    pal = brewer.pal(K, "Dark2")
    df = data.frame(x = m, y = u)
    col = as.factor(class_i)
    levels(col) = c("Dark green", "Orange", "Purple", "Pink", 
                    "Light green", "Yellow", "Brown", "Gray")[1:K]
    plot(df$x, df$y, col = as.character(col), pch = 19, main = CpG, 
         xlab = "M", ylab = "U", xlim = c(0, max(df$x) + 100), 
         ylim = c(0, max(df$y) + 100))
    return(as.character(col))
  }
  else if (K == 2) {
    df = data.frame(x = m, y = u)
    coefA = aggregate(. ~ class_i, df[c("x")], mean)$x
    coefB = aggregate(. ~ class_i, df[c("y")], mean)$y
    pos1 = which.min((coefA + coefB)/2)
    pos2 = which.max((coefA + coefB)/2)
    real_order <- c(pos1, pos2)
    class_i <- factor(class_i, levels = c(1, 2, 3, 4))
    levels(class_i) <- as.character(match(1:4, real_order))
    col_vec <- as.vector(class_i)
    out = col_vec
    col_vec[col_vec == "1"] <- scales::alpha("dodgerblue3", 0.5)
    col_vec[col_vec == "2"] <- scales::alpha("brown3", 0.5)
    out[out == "1"] <- "blue"
    out[out == "2"] <- "red"
    plot(df$x, df$y, col = col_vec, pch = 19, main = CpG, 
         xlab = "M", ylab = "U", xlim = c(0, max(df$x) + 100), 
         ylim = c(0, max(df$y) + 100))
    return(out)
  }
  else if (K == 3) {
    df = data.frame(x = m, y = u)
    coefs <- numeric()
    lm.mod <- lm(df$y[class_i == 1] ~ df$x[class_i == 1])
    coefs[1] <- unname(lm.mod$coefficients[2])
    lm.mod <- lm(df$y[class_i == 2] ~ df$x[class_i == 2])
    coefs[2] <- unname(lm.mod$coefficients[2])
    lm.mod <- lm(df$y[class_i == 3] ~ df$x[class_i == 3])
    coefs[3] <- unname(lm.mod$coefficients[2])
    class_i <- factor(class_i, levels = c(1, 2, 3))
    levels(class_i) <- as.character(match(1:3, order(coefs, 
                                                     decreasing = T)))
    col_vec <- as.vector(class_i)
    out = col_vec
    col_vec[col_vec == "1"] <- scales::alpha("dodgerblue3", 0.5)
    col_vec[col_vec == "2"] <- scales::alpha("darkmagenta", 0.5)
    col_vec[col_vec == "3"] <- scales::alpha("brown3", 0.5)
    out[out == "1"] <- "blue"
    out[out == "2"] <- "purple"
    out[out == "3"] <- "red"
    plot(df$x, df$y, col = col_vec, pch = 19, main = CpG, 
         xlab = "M", ylab = "U", xlim = c(0, max(df$x) + 100), 
         ylim = c(0, max(df$y) + 100))
    return(out)
  }
  else if (K == 4) {
    df = data.frame(x = m, y = u)
    df_norm = df/data.frame(x = rep(max(df$x), nrow(df)), 
                            y = rep(max(df$y), nrow(df)))
    df_norm$label = class_i
    coefA = aggregate(. ~ label, df_norm[c("x", "label")], 
                      mean)$x
    coefB = aggregate(. ~ label, df_norm[c("y", "label")], 
                      mean)$y
    pos1 = which.min((coefA + coefB)/2)
    coefA[pos1] = NA
    coefB[pos1] = NA
    pos4 = which.min(abs(coefA - coefB))
    coefA[pos4] = NA
    coefB[pos4] = NA
    pos2 = which.max(coefB)
    pos3 = which.max(coefA)
    real_order <- c(pos1, pos2, pos3, pos4)
    class_i <- factor(class_i, levels = c(1, 2, 3, 4))
    levels(class_i) <- as.character(match(1:4, real_order))
    col_vec <- as.vector(class_i)
    out = col_vec
    col_vec[col_vec == "1"] <- scales::alpha("black", 0.5)
    col_vec[col_vec == "2"] <- scales::alpha("dodgerblue3", 0.5)
    col_vec[col_vec == "3"] <- scales::alpha("brown3", 0.5)
    col_vec[col_vec == "4"] <- scales::alpha("darkmagenta", 0.5)
    out[out == "1"] <- "blue"
    out[out == "2"] <- "purple"
    out[out == "3"] <- "red"
    out[out == "4"] <- "red"
    plot(df$x, df$y, col = col_vec, pch = 19, main = CpG, 
         xlab = "M", ylab = "U", xlim = c(0, max(df$x) + 100), 
         ylim = c(0, max(df$y) + 100))
    return(out)
  }

}

