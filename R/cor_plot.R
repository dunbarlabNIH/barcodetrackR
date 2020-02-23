#' Correlation Plot
#'
#' Gives the pairwise correlation between each sample-sample pair in the data frame.
#'
#'@param your_SE A Summarized Experiment object.
#'@param plot_labels Vector of x axis labels. Defaults to colnames(your_SE).
#'@param method_corr Character. One of "pearson", "spearman", or "kendall". Can also use "manhattan" to compute manhattan distance instead.
#'@param your_title The title for the plot.
#'@param grid Logical. Include a grid or not in the correlation plot
#'@param label_size The size of the column labels.
#'@param plot_type Character. One of "circle", "square", "ellipse", "number", "shade", "color", or "pie".
#'@param no_negatives Logical. Whether to make negative correlations = 0.
#'@param print_tables Logical. Whether or not to print tables of p-values, confidence intervals, and R values instead of displaying a plot.
#'@param color_scale Character. One of "default", "rainbow", or "white_heat"
#'@return Plots correlation plot for your_SE.
#'@examples
#'cor_plot(your_SE = zh33, your_title = "Pearson correlation of all samples", plottype = "color")
#'@export

cor_plot = function(your_SE,
                    plot_labels = colnames(your_SE),
                    method_corr ="pearson",
                    your_title = "",
                    grid = TRUE,
                    label_size = 1,
                    plot_type = "color",
                    no_negatives = FALSE,
                    print_tables = FALSE,
                    colorscale = "default") {

  #empty table of  values
  ctab = matrix(ncol = ncol(your_data), nrow = ncol(your_data))

  #empty table for computing pvalue (not really helpful but to satisfy reviewer)
  ctab_pval = matrix(ncol = ncol(your_data), nrow =  ncol(your_data))

  #empty table for computing lower and upper bounds for 95% confidence interval
  ctab_ci_lo = matrix(nrow=ncol(your_data), ncol=ncol(your_data))
  ctab_ci_hi = matrix(nrow=ncol(your_data), ncol=ncol(your_data))

  #extracts percentages assay from your_SE
  plotting_data <- tibble::rownames_to_column(SummarizedExperiment::assays(your_SE)[["percentages"]], var = "sequence")

  #calculate the pairwise correlation of all samples with two nested loops (dual zeros would putatively artifically infalte the correlations if included, which is why this nested loop function is used)
  for (i in 1:ncol(plotting_data)) {
    for (j in 1:ncol(plotting_data)) {

      #need to run through loops of pairwise correlations to exclude those rows where dual zeros exist
      temp_df <- data.frame(plotting_data[[i]],plotting_data[[j]])
      temp_df <- temp_df[(temp_df[[1]] + temp_df[[2]]) > 0, ]
      ctab[i,j] <- cor(tempdf[[1]], tempdf[[2]], method = method_corr)
      ct_results <- cor.test(tempdf[[1]], tempdf[[2]], method = method_corr)

      #this adds p-values from the correlations and confidence intervals to a dataframe as well
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

  color_limits = c(-1,1)
  if(no_negatives){
    ctab[ctab < 0] <- 0
    color_limits <- c(0,1)
  }

  colnames(ctab) <- plot_labels
  rownames(ctab) <- plot_labels
  colnames(ctab_pval) <- plot_labels
  rownames(ctab_pval) <- plot_labels
  colnames(ctab_ci_lo) <- plot_labels
  rownames(ctab_ci_lo) <- plot_labels
  colnames(ctab_ci_hi) <- plot_labels
  rownames(ctab_ci_hi) <- plot_labels

  grid_color = ifelse(show_grid, "black", "white")

  if(print_tables){
    return(list("cortable" = ctab, "cortable_pval" = ctab_pval, "cortable_ci_hi" = ctab_ci_hi, "cortable_ci_lo" = ctab_ci_lo))
  }

  corrplot::corrplot(ctab,
                     tl.cex = label_size,
                     method = plot_type,
                     title = title(your_title, line = -1, cex.main = 2),
                     addgrid.col = grid_color,
                     tl.col="black",
                     outline = TRUE,
                     mar = c(1, 1, 3, 1),
                     col = switch(colorscale,
                                  default = NULL,
                                  rainbow = colorRampPalette(rainbow(7))(100),
                                  white_heat = colorRampPalette(col = c("black", "black", "black", "black", "black", "black","brown", "red", "orange", "yellow", "white"))(100)),
                     cl.lim = colorlimits
  )
}






