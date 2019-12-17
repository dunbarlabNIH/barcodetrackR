#' Bar or line Plot of Hematopoietic Contribution from top n clones
#'
#' Usually used for tracking a cell lineage's top clones over time.
#'
#'@param your_SE A summarized experiment. With barcode data stored in the 1st assay and metadata stored in colData.
#'@param graph_type Choice of "bar" or "line" for how to display the clonal contribution data
#'@param filter_by Name of metadata column to filter by e.g. Cell_type
#'@param filter_selection The value of the filter column to display e.g. T cells
#'@param plot_by The column of metadata that you want to be the x-axis of the plot. e.g. Timepoint
#'@param n_clones Numeric. Number of top clones from each sample that should be displayed.
#'@param linesize Numeric. Thickness of the lines. OPTIONAL
#'@param text_size Numeric. Size of text in plot. OPTIONAL
#'@return Displays a stacked area plot (made by ggplot2) of the samples' top clones.
#'@examples
#'clonal_contribution(your_data = ZG66_simple_data, graph_type = "bar",filter_by = "Cell_type", filter_selection = "B", plot_by = "Timepoint", n_clones = 20)
#'@export

clonal_contribution <- function(your_SE, graph_type = "bar", filter_by, filter_selection, plot_by, n_clones = 10, linesize = 0.5, text_size = 15){
  # Extract bc data and metadata
  your_data <- SummarizedExperiment::assay(your_SE)
  meta_data <- SummarizedExperiment::colData(your_SE)
  # Only keep data that matches filter
  your_data <- your_data[,meta_data[,filter_by] == filter_selection]
  meta_data <- meta_data[meta_data[,filter_by] == filter_selection,]
  # Name data columns based on meta data to plot by
  colnames(your_data) <- meta_data[,plot_by]
  # Take proportion and rank the barcodes
  your_data <- as.data.frame(100*prop.table(as.matrix(your_data), margin = 2))
  your_data[your_data == 0] <- NA
  your_data_ranked <- apply(-your_data, 2, rank, ties.method = "min", na.last = "keep")
  your_data$EMERGENCE <- apply(your_data_ranked, 1, function(x){return(which(x <= n_clones)[1])})
  your_data[is.na(your_data)] <- 0
  # Make a row of data with the rest of the clones
  your_data <- rbind(your_data[your_data$EMERGENCE != 0,], colSums(your_data[your_data$EMERGENCE == 0,]))
  your_data <- your_data[order(-your_data$EMERGENCE),]
  your_data$ID <- rownames(your_data)
  # Melt dataframe
  melty <- reshape2::melt(your_data, id.vars = c("EMERGENCE", "ID"), value.name = "PERCENTAGE", variable.name = "Sample")
  # melty$MONTH = as.numeric(as.character(melty$MONTH))
  melty$EMERGENCE = as.factor(melty$EMERGENCE)
  
  if (graph_type == "bar"){
      ggplot2::ggplot(melty, ggplot2::aes(x=Sample, y = PERCENTAGE, group = ID, fill = EMERGENCE))+
      ggplot2::geom_bar(stat = "identity", colour = "black",size=linesize)+
      ggplot2::labs(title = paste(filter_by,"=",filter_selection))+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))+
      # gplot2::geom_area(position = "stack", linetype = 1, size = linesize, colour = "black")+
      # ggplot2::scale_x_continuous(breaks = seq(1, months[length(months)], by = 0.5), name = "Month")+
      ggplot2::scale_y_continuous(breaks = seq(0,100,by = 10), name = "Clonal Contribution %")+
      #switch(ggplot2::plot_theme, BW = ggplot2::theme_bw(), classic = ggplot2::theme_classic(), original = ggplot2::theme_grey())+
      ggplot2::scale_fill_manual(values = c("grey", rainbow(length(unique(melty$EMERGENCE))-1,
                                                            s = 1, v = 1, start = 0, end = 0.75, alpha = 1)))+
      ggplot2::theme(text = ggplot2::element_text(size=text_size))
    # ggplot2::coord_cartesian(ylim = c(0, y_limit))
  }

  else if (graph_type == "line"){
    ggplot2::ggplot(melty, ggplot2::aes(x=Sample, y = PERCENTAGE, group = ID, fill = EMERGENCE))+
    ggplot2::geom_area(position = "stack", linetype = 1, size = linesize, colour = "black")+
    ggplot2::labs(title = paste(filter_by,"=",filter_selection))+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))+
    # ggplot2::scale_x_continuous(breaks = seq(1, months[length(months)], by = 0.5), name = "Month")+
    ggplot2::scale_y_continuous(breaks = seq(0,100,by = 10), name = "Clonal Contribution %")+
    #switch(ggplot2::plot_theme, BW = ggplot2::theme_bw(), classic = ggplot2::theme_classic(), original = ggplot2::theme_grey())+
    ggplot2::scale_fill_manual(values = c("grey", rainbow(length(unique(melty$EMERGENCE))-1,
                                       s = 1, v = 1, start = 0, end = 0.75, alpha = 1)))+
    ggplot2::theme(text = ggplot2::element_text(size=text_size))
    # ggplot2::coord_cartesian(ylim = c(0, y_limit))
  }

}
