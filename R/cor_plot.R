#' Correlation Plot
#'
#' Gives the pairwise correlation between each sample-sample pair in the data frame.
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param names Character vector. Labels for each sample.
#'@param thresh Numeric. Threshold applied to your_data and scaled for each sample. When finding correlation,
#'        only barcodes whose sum in pair of samples being compared is greater than the average of each's
#'        scaled threshold will be used for finding correlation.
#'@param your_title Title for plot.
#'@param method_corr Character. One of "pearson", "spearman", or "kendall".
#'@param labelsizes Numeric. Size of the plot labels.
#'@param plottype Character. One of "circle", "square", "ellipse", "number", "shade", "color", or "pie".
#'@param writetables Logical. Whether or not to write p-values, confidence intervals, and R values to a zip file
#'        in the current directory.
#'@param nonegatives Logical. Whether to make negative correlations = 0. Should almost always be TRUE.
#'@return Plots correlation plot of the data frame.
#'@examples
#'cor_plot(your_data = zh33, thresh = 874, your_title = "All Samples", plottype = "ellipse")
#'@export

cor_plot = function(your_data, names=colnames(your_data), thresh = 0, your_title = "", method_corr ="pearson", labelsizes = 1, plottype = "color", writetables = FALSE, no_negatives = FALSE) {

  #scales thresh to each column (sample)
  threshes <- thresh/colSums(your_data)

  #scales all the data down to a percent of its samples reads
  your_data <- as.data.frame(prop.table(as.matrix(your_data), margin = 2))
  #note - needs to be coerced into a data frame to work appropriately

  #empty table of r values
  ctab = matrix(ncol = ncol(your_data), nrow = ncol(your_data))

  #empty table for computing pvalue (not really helpful but to satisfy reviewer)
  ctab_pval = matrix(ncol = ncol(your_data), nrow =  ncol(your_data))

  #empty table for computing lower and upper bounds for 95% confidence interval
  ctab_ci_lo = matrix(nrow=ncol(your_data), ncol=ncol(your_data))
  ctab_ci_hi = matrix(nrow=ncol(your_data), ncol=ncol(your_data))

  #compare each cell type, excluding any barcodes whose sum in the two samples being compared
  #is less than the average of their two thresholds
  for (i in 1:ncol(your_data)) {
    for (j in 1:ncol(your_data)) {
      tempdf <- data.frame(your_data[[i]],your_data[[j]])
      tempdf <- tempdf[(tempdf[[1]] + tempdf[[2]]) > ((threshes[i] + threshes[j])/2), ]
      ctab[i,j] <- cor(tempdf[[1]], tempdf[[2]], method = method_corr)


      ct_results <- cor.test(tempdf[[1]], tempdf[[2]], method = method_corr)

      if(is.null(ct_results$p.value)){
        ctab_pval[i,j] <- "NA"
      } else {
        ctab_pval[i,j] <- ct_results$p.value
      }

      if(is.null(ct_results$conf.int)){
        ctab_ci_lo[i,j] = "NA"
        ctab_ci_hi[i,j] = "NA"
      } else {
        ctab_ci_lo[i,j]=ct_results$conf.int[1]
        ctab_ci_hi[i,j]=ct_results$conf.int[2]
      }

    }
  }

  ctab[is.na(ctab)] <- 0
  colorlimits = c(-1,1)
  if(no_negatives == TRUE){
    ctab[ctab < 0] <- 0
    colorlimits <- c(0,1)
  }

  colnames(ctab) <- names
  rownames(ctab) <- names
  colnames(ctab_pval) <- names
  rownames(ctab_pval) <- names
  colnames(ctab_ci_lo) <- names
  rownames(ctab_ci_lo) <- names
  colnames(ctab_ci_hi) <- names
  rownames(ctab_ci_hi) <- names

  corrplot::corrplot(ctab,tl.cex = labelsizes, method=plottype, title = title(your_title, line = -1, cex.main = 2), addgrid.col="white", tl.col="black", outline=TRUE, col=colorRampPalette(rainbow(7))(200), mar = c(1, 1, 3, 1), cl.lim = colorlimits)


  if (writetables == TRUE){
    write.table(ctab, paste(your_title,"_cortable.txt", sep = ""), quote = FALSE, sep = '\t')
    write.table(ctab_pval, paste(your_title,"_cortab_pval.txt", sep = ""), quote = FALSE, sep = '\t')
    write.table(ctab_ci_lo, paste(your_title,"_cortab95%ci_low.txt", sep = ""), quote = FALSE, sep = '\t')
    write.table(ctab_ci_hi, paste(your_title,"_cortab95%ci_hi.txt", sep = ""), quote = FALSE, sep = '\t')
  }


}
