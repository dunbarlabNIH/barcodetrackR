#' Diversity Plot
#'
#' Plots (using ggplot2) a line plot that tracks the diversity of a lineage of cells across time points.
#'
#'@param your_SE A summarized experiment. With barcode data stored in the 1st assay and metadata stored in colData.
#'@param index_type Character. One of "herfindahl-hirschman", "shannon", "simpson", or "invsimpson".
#'@param measure The measure you want to plot. Can be "diversity", "count", or "evenness"
#'@param plot_by The column of metadata that you want to be the x-axis of the plot. e.g. timepoint
#'@param group_by The column of metadata you want to group by e.g. cell _type
#'@param point_size Numeric. Size of points. Optional
#'@param line_size Numeric. Size of lines. Optional
#'@return Outputs plot of diversity tracked over time.
#'@examples
#'diversity_plot(your_data = ZG66_simple_data, index_type = "shannon",
#'                measure = "diversity", plot_by = timepoint, group_by = cell_type)
#'
#'@export
diversity_plot <- function(your_SE, index_type = "shannon", measure = "diversity", plot_by, group_by, point_size = 5, line_size = 2) {
  # Extract bc data and metadata
  your_data <- SummarizedExperiment::assay(your_SE)
  meta_data <- SummarizedExperiment::colData(your_SE)
  
  # Make plot_by an ordered factor if it's categorical
  if (class(meta_data[,plot_by]) == "factor") {
    # Order meta_data based on the order of the data not alphabetically
    meta_data[,plot_by] <- factor(meta_data[,plot_by], levels = unique(meta_data[,plot_by]))
    meta_data <- meta_data[order(meta_data[,plot_by]),]
    your_data <- your_data[,order(meta_data[,plot_by])]
  }
  
  # Calculate diversity based on specified index
  herfindahl = FALSE
  if(index_type == "herfindahl-hirschman"){
    herfindahl = TRUE
    index_type = "simpson"
  }
  div_data <- vegan::diversity(your_data, index = index_type, MARGIN = 2)
  if(herfindahl){
    div_data <- (1-your_data)
  }
  # Calculate diversity count or evenness if specified
  if (measure == "count"){
    div_data <- exp(div_data)
  }

  if (measure == "evenness"){
    div_data <- div_data / log2(colSums(your_data))
  }
  
  # Create dataframe for plotting
  plot.df <- data.frame(div_data=div_data,
                        x_data = meta_data[,plot_by],
                        group_data = meta_data[,group_by])
  
  # Correct ordering of plot_by variable
  if (class(meta_data[,plot_by]) == "numeric"){
    plot.df <- plot.df[order(plot.df$x_data),]
    # Change to factor
    plot.df$x_data <- factor(plot.df$x_data)
  } 
  
  # Create ggplot
  ggplot2::ggplot(data=plot.df, ggplot2::aes(x=x_data, y=div_data, group=group_data, colour=group_data)) +
    ggplot2::geom_line(size = line_size)+
    ggplot2::geom_point(size = point_size)+ ggplot2::labs(y = measure, x = plot_by, col = group_by)+
    ggplot2::theme_classic() + ggplot2::theme(axis.text.x = ggplot2::element_text(face="bold", size=14,color = "black"),
                            axis.text.y = ggplot2::element_text(face="bold", size=14, color = "black"))+
    ggplot2::theme(text = ggplot2::element_text(size=20))
}
