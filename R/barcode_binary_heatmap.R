#' Barcode Binary Heatmap
#'
#' Creates a binary heatmap showing the absence or presence of new clones in samples ordered from L to R in the SummarizedExperiment.
#'
#'@param your_SE A Summarized Experiment object.
#'@param plot_labels Vector of x axis labels. Defaults to colnames(your_SE).
#'@param threshold Clones with a proportion below this threshold will be set to 0.
#'@param your_title The title for the plot.
#'@param label_size The size of the column labels.
#'@param return_table  Logical. Whether or not to return table of barcode sequences with their presence or absence in each sample indicated as a 1 or 0 resepctively in the value column column.
#'@return Displays a binary heat map in the current plot window. Or if return_table is set to TRUE, returns a dataframe indicating the presence or absence of each barcode in each sample. 
#'
#'@export
#'
#'@examples
#'barcode_binary_heatmap(your_SE = wu_subset[,1:4])
#'
barcode_binary_heatmap <- function(your_SE,
                                   plot_labels = NULL,
                                   threshold = 0,
                                   your_title = NULL,
                                   label_size = 12,
                                   return_table = FALSE) {

  #get labels for heatmap
  plot_labels <- plot_labels %||% colnames(your_SE)
  if(length(plot_labels) != ncol(your_SE)){
    stop("plot_labels must be same length as number of columns being plotted")
  }

  #changing plotting_data to binary
  plotting_data <- SummarizedExperiment::assays(your_SE)[["proportions"]]
  plotting_data[plotting_data < threshold] <- 0
  plotting_data[plotting_data > 0] <- 1
  plotting_data <- plotting_data[rowSums(plotting_data) > 0,]
  barcode_order <- rownames(plotting_data)[do.call(order, plotting_data)]

  #organizing data for plotting
  plotting_data <- tibble::rownames_to_column(plotting_data, var = "sequence")
  plotting_data <- tidyr::pivot_longer(plotting_data, cols = -sequence, names_to = "sample_name", values_to = "value")
  plotting_data$sample_name <- factor(plotting_data$sample_name, levels = colnames(your_SE))
  plotting_data$value <- factor(plotting_data$value, levels = c(0,1))
  plotting_data$sequence <- factor(plotting_data$sequence, levels = barcode_order)
  x_column <- plyr::mapvalues(plotting_data$sample_name, from = colnames(your_SE), to = plot_labels)
  

  if(return_table){
    return(plotting_data)
  }
  
  plotting_data$x_labels <- factor(x_column, levels = plot_labels)

  ggplot2::ggplot(plotting_data, ggplot2::aes(x = x_labels, y = sequence))+
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
      plot.title = ggplot2::element_text(size = label_size),
      axis.text.x = ggplot2::element_text(angle=90, hjust = 1, vjust = 0.5, size = label_size),
      legend.title = ggplot2::element_text(size =  label_size),
      legend.key.width=ggplot2::unit(0.2, "cm"),
      legend.text = ggplot2::element_text(size = label_size),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1)
    )


}
