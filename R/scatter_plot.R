#' Scatter Plot
#'
#' Plots a scatter plot of two samples in the Summarized Experiment object
#'
#'@param your_SE A Summarized Experiment object of two samples.
#'@param assay The choice of assay to plot on the scatter plot
#'@param plot_labels The labels for the X and Y axis of the plot
#'@param method_corr Character. One of "pearson", "spearman", or "kendall". Can also use "manhattan" to compute manhattan distance instead.
#'@param display_corr Logical. Whether to display the computer correlation or not.
#'@param point_size The size of the points being plotted.
#'@param your_title The title for the plot.
#'@param text_size Numeric. Size of text in plot.
#'@return Plots scatter plot for the samples in your_SE.
#'
#'@importFrom rlang %||%
#'@importFrom magrittr %>%
#'
#'@export
#'
#'@examples
#'scatter_plot(your_SE = zh33[,1:2])
#"

scatter_plot = function(your_SE,
                        assay = "percentages",
                        plot_labels = colnames(your_SE),
                        method_corr ="pearson",
                        display_corr = TRUE,
                        point_size = 0.5,
                        your_title = "",
                        text_size = 12) {

  #extracts assay from your_SE
  plotting_data <- SummarizedExperiment::assays(your_SE)[[assay]]
  if(ncol(plotting_data) != 2){
    stop("your_SE must only contain 2 samples (columns)")
  }


  colnames(plotting_data) <- plot_labels

  plotting_data <- plotting_data[rowSums(plotting_data) > 0,]
  if(display_corr){
    if (method_corr == "manhattan"){
      correlation_label <- paste0(method_corr," dist: ", round(cor.test(plotting_data[[1]], plotting_data[[2]], method = method_corr)$estimate,3))
    } else {
      correlation_label <- paste0(method_corr," corr: ", round(cor.test(plotting_data[[1]], plotting_data[[2]], method = method_corr)$estimate,3))
    }
      } else {
    correlation_label <- NULL
  }

  gg_scatter <- ggplot2::ggplot(plotting_data, ggplot2::aes_(x = as.name(colnames(plotting_data)[1]),
                                                             y = as.name(colnames(plotting_data)[2]))) +
    ggplot2::geom_point(size = point_size, color = 'black')+
    ggplot2::theme(panel.background  = ggplot2::element_rect(color = "black", fill = "white"), panel.grid = ggplot2::element_blank())+
    ggplot2::xlab(plot_labels[1])+
    ggplot2::ylab(plot_labels[2])+
    ggplot2::ggtitle(paste0(your_title, "\n", correlation_label))+
    ggplot2::theme(text = ggplot2::element_text(size=text_size))

  gg_scatter

}






