#' Barcode Binary Heatmap
#'
#' Creates a binary heatmap showing the presence of new clones from L to R in the dataset.
#'
#'@param your_SE A Summarized Experiment object.
#'@param plot_labels Vector of x axis labels. Defaults to colnames(your_SE).
#'@param threshold Clones with a proportion below this threshold will be set to 0.
#'@param your_title The title for the plot.
#'@param label_size The size of the column labels.
#'@return Displays a binary heat map in the current plot window.
#'
#'@export
#'
#'@examples
#'barcode_binary_heatmap(your_SE = ZH33_SE[,1:10])
#'
barcode_binary_heatmap <- function(your_SE,
                                   plot_labels = NULL,
                                   threshold = 0,
                                   your_title = "",
                                   label_size = 12) {

  #get labels for heatmap
  plot_labels <- plot_labels %||% colnames(your_SE)
  if(length(plot_labels) != ncol(your_SE)){
    stop("plot_labels must be same length as number of columns being plotted")
  }

  #changing plotting_data to binary
  plotting_data <- SummarizedExperiment::assays(your_SE)[["percentages"]]
  plotting_data[plotting_data < threshold] <- 0
  plotting_data[plotting_data > 0] <- 1
  plotting_data <- plotting_data[rowSums(plotting_data) > 0,]
  barcode_order <- rownames(plotting_data)[do.call(order, plotting_data)]

  #organizing data for plotting
  plotting_data <- tibble::rownames_to_column(plotting_data, var = "sequence")
  plotting_data <- tidyr::pivot_longer(plotting_data, cols = -sequence, names_to = "sample_name", values_to = "value")
  plotting_data$sample_name <- factor(plotting_data$sample_name, levels = plot_labels)
  plotting_data$value <- factor(plotting_data$value, levels = c(0,1))
  plotting_data$sequence <- factor(plotting_data$sequence, levels = barcode_order)


  ggplot2::ggplot(plotting_data, ggplot2::aes(x = sample_name, y = sequence))+
    ggplot2::geom_tile(ggplot2::aes(fill = value))+
    ggplot2::scale_fill_manual(paste0("Detection"),
                               values = c("0" = "white", "1" = "#4575B4"),
                               expand = c(0,0),
                               labels = c("not detected", "detected"))+
    ggplot2::scale_y_discrete(labels = NULL, breaks = NULL, expand = c(0,0))+
    ggplot2::scale_x_discrete(expand = c(0,0), labels = plot_labels)+
    ggplot2::ylab(NULL)+
    ggplot2::xlab(NULL)+
    ggplot2::ggtitle(your_title)+
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 20),
      axis.text.x = ggplot2::element_text(angle=90, hjust = 1, vjust = 0.5, size = label_size),
      legend.title = ggplot2::element_text(size =  label_size, face = "bold"),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1)
    )


}
