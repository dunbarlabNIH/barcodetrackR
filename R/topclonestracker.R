#' Top Clones Tracker
#'
#' Given a chosen sample, will track the top clones from that samples over the other samples included in the dataset.
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param your_data A data frame of one column, with rownames being the barcodes. Used to get the top clones.
#'@param months Numeric vector. Vector corresponding to the months of each sample.
#'@param n_clones Numeric. Number of top clones from each month that should be displayed.
#'@param linesize Numeric. Thickness of the lines.
#'@param text_size Numeric. Size of text in plot.
#'@param y_limit Numeric. Limit of y axis in plot.
#'@param plot_theme Character. One of "BW", "classic", or "original".
#'@param your_title Your title for the plot.
#'@return Displays a stacked area plot (made by ggplot2) of the samples' top clones.
#'@examples
#'topclonestracker(your_data = zh33_Tlineage, top_clones_choice = zh33T1m, months = c(1,2,3,4,5,6.5, 9.5, 12), y_limit = 50)
#'@export

topclonestracker <- function(your_data, top_clones_choice, months = c(1:ncol(your_data)), n_clones = 10,
                             linesize = 0.5, text_size = 15, y_limit = 100, plot_theme = "classic",
                             your_title = ""){
  top_barcodes <- rownames(top_clones_choice[order(-top_clones_choice),,drop = FALSE])[1:n_clones]
  colnames(your_data) <- months
  your_data <- as.data.frame(100*prop.table(as.matrix(your_data), margin = 2))
  your_data <- rbind(your_data[rownames(your_data) %in% top_barcodes,],
                     colSums(your_data[!(rownames(your_data) %in% top_barcodes),,drop = FALSE]))
  rownames(your_data)[nrow(your_data)] <- "OTHERS"
  your_data <- reshape2::melt(as.matrix(your_data))
  colnames(your_data) <- c("BARCODE", "MONTH", "PROP")
  ggplot2::ggplot(your_data, ggplot2::aes(x=MONTH, y = PROP, group = BARCODE, fill = BARCODE))+
    ggplot2::geom_area(position = "stack", linetype = 1, size = linesize, colour = "black")+
    ggplot2::scale_x_continuous(breaks = seq(1, months[length(months)], by = 0.5), name = "Month")+
    ggplot2::scale_y_continuous(breaks = seq(0,100,by = 10), name = "Hematopoietic Contribution %")+
    switch(plot_theme, BW = ggplot2::theme_bw(), classic = ggplot2::theme_classic(), original = ggplot2::theme_grey())+
    ggplot2::scale_fill_manual(guide = FALSE, values = c(rainbow(length(unique(your_data$BARCODE))-1,
                                                                 s = 1, v = 1, start = 0, end = 0.75, alpha = 1), "grey"))+
    ggplot2::theme(text = ggplot2::element_text(size=text_size))+
    ggplot2::ggtitle(your_title)+
    ggplot2::coord_cartesian(ylim = c(0, y_limit))

}
