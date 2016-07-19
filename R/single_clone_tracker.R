#' Single Clone Tracker
#'
#' Given a data frame organized from first timepoint to latest, plots top N clones over time.
#'@param your_data Your data frame.
#'@param n_clones Number of clones to choose.
#'@param y_max Maximum value of plot.
#'@param line_size Size of the lines.
#'@param top_value Max point at which to set the gradient to red.
#'@param text_size Text size for the whole plot.
#'@return Line plot for a chosen lineage over time.
#'@examples
#'single_clone_tracker(datum[,c("T 1m", "T 3m" "T 5m")], y_max = 0.1, top_value = 0.05)
#'@export
single_clone_tracker <- function(your_data, n_clones = 100, y_max = 1, line_size = 2, top_value = .5, text_size = 15) {
  your_data <- as.data.frame(prop.table(as.matrix(your_data), margin = 2))
  your_data <- your_data[barcodetrackR::gettopindices(your_data, top = n_clones),]
  colnames(your_data) <- 1:ncol(your_data)
  your_data$TOTAL <- rowSums(your_data)
  your_data$TOTAL <- apply(your_data,1,max)
  your_data <- your_data[order(your_data$TOTAL),]
  your_data$GROUP <- as.factor(1:nrow(your_data))
  your_data <- reshape2::melt(your_data, id.vars = c("TOTAL", "GROUP"))
  ggplot2::ggplot(your_data, ggplot2::aes(x = variable, y = value, color = TOTAL, group = GROUP))+
    ggplot2::geom_line(size = line_size, lineend = 'round')+
    ggplot2::scale_y_continuous(labels = function(x){paste(round(x*100, 2), "%")})+
    ggplot2::scale_color_gradientn(colors = colorRampPalette(colors = c("grey", "red", "red"))(100), name = "Abundance\n", values = c(0,top_value,1), limits = c(0,1), labels = function(x){paste(x*100, "%")})+
    ggplot2::coord_cartesian(ylim = c(0,y_max))+
    ggplot2::ylab("\nAbundance\n")+
    ggplot2::xlab("\nTimepoint\n")+
    ggplot2::theme(text = ggplot2::element_text(size = text_size),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(colour = "grey"),
                   panel.background = ggplot2::element_rect(fill = "white", colour = "black"),
                   panel.margin = ggplot2::unit(2,"lines"),
                   axis.text.x = ggplot2::element_text(vjust = 0.5, hjust = 1, angle = 90),
                   axis.ticks.length = ggplot2::unit(.75, "cm"),
                   axis.ticks = ggplot2::element_line(size = 1))


}
