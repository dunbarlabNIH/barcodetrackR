#'@importFrom rlang %||%
#'@importFrom magrittr %>%
#'
#'@title Bias histogram
#'
#'@description Given a summarized experiment, gives histogram of log biases for 2 cell types
#'
#'@param your_SE Your SummarizedExperiment of barcode data and associated metadata
#'@param split_bias_on The column in `colData(your_SE)` from which `bias_1` and `bias_2` will be chosen
#'@param bias_1 The factor you wish to plot on the right side of the plots.
#'@param bias_2 The factor you wish to plot on the left side of the plots.
#'@param split_bias_over The column in `colData(your_SE)` that you wish to observe the bias split on.
#'@param bias_over Choice(s) from the column designated in `split_bias_over` that will be used for plotting. Defaults to all.
#'@param breaks Numeric. The breaks specified for bins on the x-axis (how biased the clones are towards one factor or the other).
#'@param text_size The size of the text in the plot.
#'@param linesize The linewidth of the stacked bars which represent individual barcodes
#'@param ncols Numeric. Number of columns to plot on using plot_grid from cowplot.
#'@param same_all_y Logical. Whether or not to plot all plots on the same y axis limits.
#'@return Histogram of log bias for two factors over another set of factors.
#'
#'@importFrom rlang %||%
#'
#'@examples
#'bias_histogram(your_SE = SE, split_bias_on = "Lineage", bias_1 = "B", bias_2 = "T", split_bias_over = "Month", bias_over = c(1,2,4.5,12))
#'@export
bias_histogram <- function(your_SE, split_bias_on, bias_1, bias_2, split_bias_over, bias_over = NULL, breaks = c(10,2,1,0.5), text_size = 10, linesize = .4, ncols = 1, scale_all_y = TRUE) {

  breaks_labels <- breaks
  breaks <- c(-Inf, sort(unique(c(-breaks, breaks))), Inf)

  # Some basic error checking before running the function
  if(length(breaks) != length(unique(breaks))){
    stop("breaks must be unique")
  }
  coldata_names <- colnames(SummarizedExperiment::colData(your_SE))
  if(any(! c(split_bias_on, split_bias_over) %in% coldata_names)){
    stop("split_bias_on and split_bias_over must both match a column name in colData(your_SE)")
  }
  if(any(! c(bias_1, bias_2) %in% unique(SummarizedExperiment::colData(your_SE)[[split_bias_on]]))){
    stop("bias_1 and bias_2 must both be elements in the colData column specified with split_bias_on")
  }
  if(is.numeric(SummarizedExperiment::colData(your_SE)[[split_bias_over]])){
    bias_over <- bias_over %||% sort(unique(SummarizedExperiment::colData(your_SE)[[split_bias_over]]))
  } else {
    bias_over <- bias_over %||% levels(SummarizedExperiment::colData(your_SE)[[split_bias_over]])
  }

  plot_list <- lapply(1:length(bias_over), function(i){
    temp_subset <- your_SE[,(your_SE[[split_bias_over]] %in% bias_over[i]) & (your_SE[[split_bias_on]] %in% c(bias_1, bias_2))]
    if(ncol(temp_subset) != 2){
      stop("Error in subsetting")
    }
    your_data <- cbind(SummarizedExperiment::assays(temp_subset[, temp_subset[[split_bias_on]] == bias_1])[["percentages"]],
                       SummarizedExperiment::assays(temp_subset[, temp_subset[[split_bias_on]] == bias_2])[["percentages"]])
    colnames(your_data) <- c("bias_1", "bias_2")
    your_data <- your_data[rowSums(your_data) > 0,] %>%
      tibble::rownames_to_column(var = "barcode") %>%
      dplyr::mutate(added_percentages = bias_1 + bias_2, bias = bias_1/bias_2) %>%
      dplyr::mutate(log2_bias = log2(bias)) %>%
      dplyr::mutate(log2_bias_cuts = cut(log2_bias,
                                         breaks = breaks,
                                         include.lowest = TRUE))
    g <- ggplot2::ggplot(your_data[order(your_data$added_percentages),],
                         ggplot2::aes(x = log2_bias_cuts, y = added_percentages))+
      ggplot2::geom_bar(stat = "identity", fill = "white", size = linesize, color = "black")+
      ggplot2::scale_x_discrete(name = paste0("log bias: log2(", bias_1, "/", bias_2, ")"), drop = FALSE)+
      ggplot2::scale_y_continuous(name = "Added Proportions", labels = function(i)(paste0(i*100, "%")))+
      ggplot2::theme(text = ggplot2::element_text(size = text_size),
                     panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_line(colour = "grey"),
                     panel.background = ggplot2::element_rect(fill = "white", colour = "black"),
                     panel.spacing = ggplot2::unit(2, "lines"),
                     axis.text.x = ggplot2::element_text(vjust = 0.5, hjust = 1, angle = 90))+
      ggplot2::labs(title = paste0(split_bias_over,": ", bias_over[i]))+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  })
  if(scale_all_y){
    plot_max <- max(unlist(lapply(plot_list, function(x){ggplot2::ggplot_build(x)$layout$panel_params[[1]]$y.range[2]})))
    plot_list <- lapply(plot_list, function(x){x + ggplot2::coord_cartesian(ylim = c(0, plot_max))})
  }
  cowplot::plot_grid(plotlist = plot_list, ncol = ncols)
}


