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
  myPalette <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", 
                 "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", 
                 "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", 
                 "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
                 "#8A7C64", "#599861")
  colnames(your_data) <- 1:ncol(your_data)
  your_data <- as.data.frame(prop.table(as.matrix(your_data), margin = 2))
  your_data <- your_data[barcodetrackR::gettopindices(your_data, top = n_clones),]
  your_data <- your_data[order(-apply(your_data,1,max)),]
  your_data$BARC <- as.factor(rev(1:nrow(your_data)))
  your_data$COLOOR <- c(RColorBrewer::brewer.pal(8, 'Set1'), "cyan", "black", rep('grey', nrow(your_data)-10))
  your_data <- reshape2::melt(your_data, id.vars = c("BARC", "COLOOR"))
  print(head(your_data, n = 100))
  print(str(your_data))
  ggplot2::ggplot(your_data, ggplot2::aes(x = variable, y = value))+
    ggplot2::geom_line(ggplot2::aes(group = BARC, color = COLOOR), size = line_size, lineend = 'round')+
    ggplot2::scale_y_continuous(labels = function(x){paste(round(x*100, 2), "%")})+
    ggplot2::scale_color_identity()+
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
