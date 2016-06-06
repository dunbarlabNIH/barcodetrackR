#' Hematopoietic Contribution
#'
#' Usually used for tracking a cell lineage's top clones over time. Best to
#' put samples in order from left to right (i.e. "1mT", "2mT", "3mT", etc.)
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param months Numeric vector. Vector corresponding to the months of each sample.
#'@param n_clones Numeric. Number of top clones from each month that should be displayed.
#'@param linesize Numeric. Thickness of the lines.
#'@param text_size Numeric. Size of text in plot.
#'@param y_limit Numeric. Limit of y axis in plot.
#'@param plot_theme Character. One of "BW", "classic", or "original".
#'@return Displays a stacked area plot (made by ggplot2) of the samples' top clones.
#'@examples
#'hematoper(your_data = zh33_Tlineage, months = c(1,2,3,4,5,6.5, 9.5, 12), y_limit = 50)
#'@export

hematoper <- function(your_data, months, n_clones = 10, linesize = 0.5, text_size = 15, y_limit = 100, plot_theme = "classic"){
  colnames(your_data) <- months
  your_data <- as.data.frame(100*prop.table(as.matrix(your_data), margin = 2))
  your_data[your_data == 0] <- NA
  your_data_ranked <- apply(-your_data, 2, rank, ties.method = "min", na.last = "keep")
  your_data$EMERGENCE <- apply(your_data_ranked, 1, function(x){return(which(x <= n_clones)[1])})
  your_data[is.na(your_data)] <- 0
  your_data <- rbind(your_data[your_data$EMERGENCE != 0,], colSums(your_data[your_data$EMERGENCE == 0,]))
  your_data <- your_data[order(-your_data$EMERGENCE),]
  your_data$ID <- rownames(your_data)
  melty <- reshape2::melt(your_data, id.vars = c("EMERGENCE", "ID"), value.name = "PERCENTAGE", variable.name = "MONTH")
  melty$MONTH = as.numeric(as.character(melty$MONTH))
  melty$EMERGENCE = as.factor(melty$EMERGENCE)
  ggplot2::ggplot(melty, ggplot2::aes(x=MONTH, y = PERCENTAGE, group = ID, fill = EMERGENCE))+
    ggplot2::geom_area(position = "stack", linetype = 1, size = linesize, colour = "black")+
    ggplot2::scale_x_continuous(breaks = seq(1, months[length(months)], by = 0.5), name = "Month")+
    ggplot2::scale_y_continuous(breaks = seq(0,100,by = 10), name = "Hematopoietic Contribution %")+
    switch(plot_theme, BW = ggplot2::theme_bw(), classic = ggplot2::theme_classic(), original = ggplot2::theme_grey())+
    ggplot2::scale_fill_manual(values = c("grey", rainbow(length(unique(melty$EMERGENCE))-1,
                                       s = 1, v = 1, start = 0, end = 0.75, alpha = 1)))+
    ggplot2::theme(text = ggplot2::element_text(size=text_size))+
    ggplot2::coord_cartesian(ylim = c(0, y_limit))

}
