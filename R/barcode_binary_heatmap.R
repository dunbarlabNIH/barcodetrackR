#' Barcode Binary Heatmap
#'
#' Creates a binary heatmap showing the presence of new clones from L to R in the dataset.
#'
#'@param your_SE A Summarized Experiment object.
#'@param plot_labels Vector of x axis labels. Defaults to colnames(your_SE).
#'@param your_threshold A chosen threshold (set as a proportion, i.e. 0.01 is 1%) that decides whether a feature entry is set to 0. Defaults to 0.
# #'@param your_title The title for the plot.
#'@param label_size The size of the column labels.
#'@return Displays a binary heat map in the current plot window.
#'@examples
#'barcode_binary_heatmap(your_SE = ZH33_SE[,1:10])
#'@export


barcode_binary_heatmap <- function(your_SE,
                                   plot_labels = NULL,
                                   your_threshold = 0,
                                  # your_title = "",
                                   label_size = 1) {

  #get labels for heatmap
  plot_labels <- plot_labels %||% colnames(your_SE)
  if(length(plot_labels) != ncol(your_SE)){
    stop("plot_labels must be same length as number of columns being plotted")
  }

  #changing plotting_data to binary
  plotting_data <- SummarizedExperiment::assays(your_SE)[["percentages"]]
  plotting_data[plotting_data < your_threshold] <- 0
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
    ggplot2::scale_fill_manual(values = c("1" = "blue", "0" = "white"), expand = c(0,0), labels = c("Absent", "Present"))+
    ggplot2::scale_y_discrete(labels = NULL, breaks = NULL, expand = c(0,0))+
    ggplot2::scale_x_discrete(expand = c(0,0), labels = plot_labels)+
    ggplot2::ylab(NULL)+
    ggplot2::xlab(NULL)+
    ggplot2::ggtitle(paste0("\n", your_title, "\n"))+
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 20),
      axis.text.x = ggplot2::element_text(angle=90, hjust = 1, vjust = 0.5, size = label_size),
      legend.text = ggplot2::element_text(size =  15, face = 'bold'),
      legend.title = ggplot2::element_text(size =  15),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1)
    )


}
