#' Computes bimodality per CpG across samples in a coefficient of variation matrix
#' @description a
#' @details a
#' @param CV A coefficient of variation matrix (probes as rows, samples as columns)
#' @return A vector of coefficient of bimodality per CpG, \code{BC}
#' @examples
#' compute_BC_CV(CV)
Visualize_cometh = function(annotation, CpG, distance, L_bound, R_bound, beta_mat, probe2exclude = NULL, max_y = 11)
{
  # CpG = 'cg01026744'
  # R_bound = 4; L_bound = 4
  # distance = 50
  annotation = annotation[order(annotation$chr, annotation$pos),]
  sub_annot = annotation[CpG,]
  tmp = annotation[annotation$chr == sub_annot$chr,]
  tmp = tmp[(tmp$pos > sub_annot$pos - distance) & (tmp$pos < sub_annot$pos + distance),]

  if(nrow(tmp) == 1)
  {
    stop('Increase distance to make plot. No other CpG could be found in the proposed window')
  }

  pos_vec = match(rownames(tmp), rownames(annotation))
  if(pos_vec[1] - L_bound <= 0)
  {
    stop('Shorten L_bound; There are not so many CpGs on the left')
  }

  if(annotation[pos_vec[length(pos_vec)],'chr'] != annotation[pos_vec[length(pos_vec)] + R_bound,'chr'])
  {
    stop('Shorten R_bound; You are already on the next chromosome')
  }

  sub_annot = annotation[(pos_vec[1]-L_bound):(pos_vec[length(pos_vec)] + R_bound),]
  #beta_mat = M[rownames(sub_annot),]/(M[rownames(sub_annot),]+U[rownames(sub_annot),]+100)
  beta_mat = beta_mat[rownames(sub_annot),]
  # beta_mat = Beta[rownames(sub_annot),]

  dat = as.matrix(cor(t(beta_mat)))^2
  rownames(dat) <- sub_annot$pos
  colnames(dat) <- sub_annot$pos



  par(mar=c(8.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  phic = plotHic(dat, sub_annot$chr, as.integer(rownames(dat)[1]), as.integer(rownames(dat)[nrow(dat)]),
                 max_y, palette = function(x) colorRampPalette(c("aliceblue", "darkblue"))(x),
                 flip = FALSE)
  #arrows(x0 = start, y0 = -0.1, x1 = end, y1 = -0.1, code = 0, xpd = T)
  addlegend(phic[[1]], palette=phic[[2]], title="", side="right", bottominset=0.4, topinset=0,
            xoffset=0.025, labelside="left", width=0.025, title.offset=0.07)

  # labelgenome(sub_annot$chr, as.integer(rownames(dat)[1]), as.integer(rownames(dat)[nrow(dat)]),
  #             side=1, scipen=20, n=0, scale="Mb", edgeblankfraction=0.30, line=.18, chromline=.5, scaleline=0.5)
  #
  start = sub_annot$pos[1]
  end = sub_annot$pos[nrow(dat)]
  jump = (end-start)/nrow(dat)
  axis(side = 1, at = seq(start+jump/2, end-jump/2, jump), labels = rep('', nrow(dat)), srt = 90)
  col_vec = rep('black', times = nrow(dat))
  if(!is.null(probe2exclude))
  {
    col_vec[rownames(sub_annot) %in% probe2exclude] = 'magenta4'
  }
  col_vec[rownames(sub_annot) == CpG] = 'red3'



  text(x = seq(start+jump/2, end-jump/2, jump), y = rep(-1, nrow(dat)), xpd = T, labels = rownames(sub_annot),
       srt = 90, cex = 0.6, col = col_vec)

  #axis(side = 1, at = seq(start+jump, end-jump+1, jump), labels = rep('', nrow(dat)-1), srt = 90, pos = c(-1.7), col = 'cornsilk4')

  deltaD = sapply(1:(nrow(dat)-1), function(x) sub_annot$pos[x+1]-sub_annot$pos[x])
  text(x = seq(start+jump, end-jump+1, jump), y = rep(-0.2, nrow(dat)), xpd = T, labels = deltaD,
       srt = 0, cex = 0.6, col = 'cornsilk4')
  par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  text(x = end, y = -0.2, labels = expression(paste(Delta, delta)), xpd = T)
  text(x = (start + end)/2, y = -2.1, labels = sub_annot$chr, xpd = T)

  cor_vec = sapply(1:nrow(beta_mat), function(x) cor(beta_mat[x, seq(1, ncol(beta_mat), 2)], beta_mat[x, seq(2, ncol(beta_mat), 2)],
                                                     method = "pearson")^2)
  names(cor_vec) = rownames(beta_mat)

  return(rownames(beta_mat))
}
