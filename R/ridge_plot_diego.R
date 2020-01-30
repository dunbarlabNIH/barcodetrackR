#'@importFrom rlang %||%
#'@importFrom magrittr %>%
#'
#'@title  Ridge plot
#'
#'@description Given a summarized experiment, gives ridge plots showing percent total contribution to both lineages.
#'
#'@param split_bias_on The column in `colData(your_SE)` from which `bias_1` and `bias_2` will be chosen
#'@param bias_1 The factor you wish to plot on the left side of the plots.
#'@param bias_2 The factor you wish to plot on the right side of the plots.
#'@param split_bias_over The column in `colData(your_SE)` that you wish to observe the bias split on.
#'@param bias_over Choice(s) from the column designated in `split_bias_over` that will be used for plotting. Defaults to all.
#'@param weighted Logical. Whether ot not to use weights when running the kernel density estimator.
#'@param text_size The size of the text in the plot.
#'@param linesize The linewidth of the stacked bars which represent individual barcodes
#'@param ncols Numeric. Number of columns to plot on using plot_grid from cowplot.
#'@return Ridge plot of log bias for two factors over another set of factors.
#'@examples
#'ridge_plot(your_SE = SE, split_bias_on = "Lineage", bias_1 = "B", bias_2 = "T", split_bias_over = "Month", bias_over = c(1,2,4.5,12))
#'@export
ridge_plot <- function(your_SE, split_bias_on, bias_1, bias_2, split_bias_over, bias_over = NULL, breaks = c(10,2,1,0.5), text_size = 10, linesize = .4, ncols = 1, weighted = TRUE, scale = 1) {

  # Some basic error checking before running the function
  coldata_names <- colnames(SummarizedExperiment::colData(your_SE))
  if(any(! c(split_bias_on, split_bias_over) %in% coldata_names)){
    stop("split_bias_on and split_bias_over must both match a column name in colData(your_SE)")
  }
  if(any(! c(bias_1, bias_2) %in% unique(SummarizedExperiment::colData(your_SE)[[split_bias_on]]))){
    stop("bias_1 and bias_2 must both be elements in the colData column specified with split_bias_on")
  }
  if(is.numeric(SummarizedExperiment::colData(your_SE)[[split_bias_over]])){
    bias_over <- bias_over %||% sort(unique(SummarizedExperiment::colData(SE)[[split_bias_over]]))
  } else {
    bias_over <- bias_over %||% levels(SummarizedExperiment::colData(your_SE)[[split_bias_over]])
  }


  subset_SE <- your_SE[,(your_SE[[split_bias_over]] %in% bias_over) & (your_SE[[split_bias_on]] %in% c(bias_1, bias_2))]
  your_data <- lapply(1:length(bias_over), function(i){
    temp_subset <- your_SE[,(your_SE[[split_bias_over]] %in% bias_over[i])]
    temp_data <- cbind(SummarizedExperiment::assays(temp_subset[, temp_subset[[split_bias_on]] == bias_1])[["normalized"]],
                       SummarizedExperiment::assays(temp_subset[, temp_subset[[split_bias_on]] == bias_2])[["normalized"]])
    colnames(temp_data) <- c("bias_1", "bias_2")
    temp_data <- temp_data[rowSums(temp_data) > 0,] %>%
      tibble::rownames_to_column(var = "barcode") %>%
      dplyr::mutate(added_normalized_values = bias_1 + bias_2, bias = ((bias_1+1)/(bias_2+1))) %>%
      dplyr::mutate(log2_bias = log2(bias)) %>%
      dplyr::mutate(!! split_bias_over := bias_over[i])
  }) %>% do.call(rbind, .)
  return(your_data)



  plot_list <- lapply(1:length(bias_over), function(i){
    temp_subset <- your_SE[,(your_SE[[split_bias_over]] %in% bias_over[i]) & (your_SE[[split_bias_on]] %in% c(bias_1, bias_2))]
    if(ncol(temp_subset) != 2){
      stop("Error in subsetting")
    }
    your_data <- cbind(SummarizedExperiment::assays(temp_subset[, temp_subset[[split_bias_on]] == bias_1])[["normalized"]],
                       SummarizedExperiment::assays(temp_subset[, temp_subset[[split_bias_on]] == bias_2])[["normalized"]])
    colnames(your_data) <- c("bias_1", "bias_2")
    your_data <- your_data[rowSums(your_data) > 0,] %>%
      tibble::rownames_to_column(var = "barcode") %>%
      dplyr::mutate(bias_1 = bias_1+1, bias_2 = bias_2+1) %>%
      dplyr::mutate(added_percentages = bias_1 + bias_2, bias = bias_1/bias_2) %>%
      dplyr::mutate(log2_bias = log2(bias))
    your_data_density <-

    g <- ggplot2::ggplot(your_data[order(your_data$added_percentages),],
                         ggplot2::aes(x = log2_bias_cuts, y = added_percentages))+
      ggplot2::geom_bar(stat = "identity", fill = "white", size = linesize, color = "black")+
      ggplot2::scale_x_discrete(name = paste0("log2 bias: log2(", bias_1, "/", bias_2, ")"), drop = FALSE)+
      ggplot2::scale_y_continuous(name = "Added Proportions", labels = function(i)(paste0(i*100, "%")))+
      ggplot2::theme(text = ggplot2::element_text(size = text_size),
                     panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_line(colour = "grey"),
                     panel.background = ggplot2::element_rect(fill = "white", colour = "black"),
                     panel.spacing = ggplot2::unit(2, "lines"),
                     axis.text.x = ggplot2::element_text(vjust = 0.5, hjust = 1, angle = 90))+
      ggplot2::labs(title = bias_over[i])+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  })




  # Weighted ridge plot
  if (weighted == T){
    as.data.frame(your_data) %>% group_by(TP) %>%
      do(ggplot2:::compute_density(.$log_bias, .$added_prop)) %>%
      dplyr::rename(log_bias = x) -> your_data_densities

    p <- ggplot2::ggplot(your_data_densities, ggplot2::aes(x = log_bias, y = TP, fill = TP, group = TP, height = density)) +
      ggridges::geom_density_ridges(stat="identity") +
      ggplot2::theme(text = ggplot2::element_text(size = text_size),
                     panel.grid.major.x = ggplot2::element_line(colour = "grey"),
                     panel.grid.major.y = ggplot2::element_blank(),
                     panel.background = ggplot2::element_rect(fill = "white"),
                     axis.text.x = ggplot2::element_text(vjust = 0.5, hjust = 1, angle = 90),
                     legend.position = "none") + ggplot2::scale_y_discrete(name = plot_by) +
                     ggplot2::annotate("text",x = Inf, y = 0.6, label = cell_1, size = 8,hjust = 5)+
                     ggplot2::annotate("text",x = -Inf, y = 0.6, label = cell_2, size = 8, hjust = -5)
  } else {

  # Normal ridge plot
  p <- ggplot2::ggplot(your_data[order(your_data$added_prop),], ggplot2::aes(x = log_bias, y = TP, fill = TP, group = TP)) +
    ggridges::geom_density_ridges(scale = scale) +
    ggplot2::theme(text = ggplot2::element_text(size = text_size),
                panel.grid.major.x = ggplot2::element_line(colour = "grey"),
                panel.grid.major.y = ggplot2::element_blank(),
                panel.background = ggplot2::element_rect(fill = "white"),
                axis.text.x = ggplot2::element_text(vjust = 0.5, hjust = 1, angle = 90),
                legend.position = "none") + ggplot2::scale_y_discrete(name = plot_by) +
                ggplot2::coord_cartesian(clip = "off")+
                ggplot2::annotate("text",x = Inf, y = -Inf, label = cell_1, size = 8,hjust = 5, vjust = -1)+
                ggplot2::annotate("text",x = -Inf, y = -Inf, label = cell_2, size = 8, hjust = -5, vjust = -1)
                # The hjust above needs to be able to scale to different x limits
  }

p

}
