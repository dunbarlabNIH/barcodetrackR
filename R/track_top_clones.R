#' Bar or line Plot of Hematopoietic Contribution from top n clones
#'
#' Usually used for tracking a cell lineage's top clones over time.
#'
#'@param your_SE A summarized experiment. With barcode data stored in the 1st assay and metadata stored in colData.
#'@param filter_by Name of metadata column to filter by e.g. Cell_type
#'@param filter_selection The value of the filter column to display e.g. T cells
#'@param plot_by The column of metadata that you want to be the x-axis of the plot. e.g. Timepoint
#'@param chosen_sample The value of the plot_by argument corresponding to the sample whose clones you want to track e.g. 82hm (82 half month)
#'@param n_clones Numeric. Number of top clones from each sample that should be displayed.
#'Note: For more than 12 clones, the colors will repeat
#'@param linesize Numeric. Thickness of the lines. OPTIONAL
#'@param text_size Numeric. Size of text in plot. OPTIONAL
#'@return Displays a stacked area plot (made by ggplot2) of the samples' top clones.
#'@examples
#'clonal_contribution(your_data = ZG66_simple_data, graph_type = "bar",filter_by = "Cell_type", filter_selection = "B", plot_by = "Timepoint", n_clones = 20)
#'@export

track_top_clones <- function(your_SE, filter_by, filter_selection, plot_by, chosen_sample, n_clones = 10, linesize = 0.5, text_size = 15){
  # Extract bc data and metadata
  your_data <- SummarizedExperiment::assay(your_SE)
  meta_data <- SummarizedExperiment::colData(your_SE)
  # Only keep data that matches filter
  your_data <- your_data[,meta_data[,filter_by] == filter_selection]
  meta_data <- meta_data[meta_data[,filter_by] == filter_selection,]
  # Order data columns based on meta data to plot by
  if (class(meta_data[,plot_by]) == "numeric"){
    meta_ordered <- meta_data[order(meta_data[,plot_by]),]
    # Change to factor
    meta_ordered[,plot_by] <- factor(meta_ordered[,plot_by])
    your_data <- your_data[,order(meta_data[,plot_by])]
  } else {
    # Order meta_data based on the order of the data not alphabetically
    meta_data[,plot_by] <- factor(meta_data[,plot_by], levels = unique(meta_data[,plot_by]))
    meta_ordered <- meta_data[order(meta_data[,plot_by]),]
    your_data <- your_data[,order(meta_data[,plot_by])]
  }
  # One basic error checking step
  if (length(which(meta_ordered[,plot_by] == chosen_sample)) == 0){
    stop("Chosen sample must be found in plot_by")
  }
  # Name data columns based on meta data to plot by
  colnames(your_data) <- meta_ordered[,plot_by]
  # Take proportion
  your_data <- as.data.frame(100*prop.table(as.matrix(your_data), margin = 2))
  your_data[your_data == 0] <- NA
  # Rank based on chosen sample
  ref <- which(meta_ordered[,plot_by] == chosen_sample)
  your_data_ordered <- your_data[order(-your_data[,ref]),]
  # return(ref)
  # return(your_data_ordered)
  #your_data_ranked <- apply(-your_data, 2, rank, ties.method = "min", na.last = "keep")
  #return(your_data_ranked)
  #your_data$EMERGENCE <- apply(your_data_ranked, 1, function(x){return(which(x <= n_clones)[1])})
  #return(your_data)
  your_data_ordered[is.na(your_data_ordered)] <- 0
  # # Make a row of data with the rest of the clones
  # # your_data <- rbind(your_data[your_data$EMERGENCE != 0,], colSums(your_data[your_data$EMERGENCE == 0,]))
  # your_data <- your_data[order(-your_data$EMERGENCE),]
  # Take only the top n clones from the chosen sample
  plot_data <- your_data_ordered[1:n_clones,]
  #return(plot_data)
  plot_data$ID <- rownames(plot_data)
  # Melt dataframe
  melty <- reshape2::melt(plot_data, id.vars = c("ID"), value.name = "PERCENTAGE", variable.name = "Sample")
  #return(melty)
  # melty$MONTH = as.numeric(as.character(melty$MONTH))
  # melty$EMERGENCE = as.factor(melty$EMERGENCE)

  # Make plot
  ggplot2::ggplot(melty, ggplot2::aes(x=Sample, y = PERCENTAGE, group = ID, fill = ID))+
    ggplot2::geom_area(position = "stack", linetype = 1, size = linesize, colour = "black")+
    ggplot2::theme_classic() +
    ggplot2::labs(title = paste(filter_by,"=",filter_selection))+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))+
    # ggplot2::scale_x_continuous(breaks = seq(1, months[length(months)], by = 0.5), name = "Month")+
    ggplot2::scale_y_continuous(name = "Clonal Contribution %")+
    #switch(ggplot2::plot_theme, BW = ggplot2::theme_bw(), classic = ggplot2::theme_classic(), original = ggplot2::theme_grey())+
    ggplot2::scale_fill_manual(values = rep(RColorBrewer::brewer.pal(12,"Paired"),times = ceiling(n_clones/12))) +
    ggplot2::theme(text = ggplot2::element_text(size=text_size))+
    ggplot2::theme(legend.position="none")
  # ggplot2::coord_cartesian(ylim = c(0, y_limit))
}
