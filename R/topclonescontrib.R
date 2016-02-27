#' Top Clones Bar Chart
#'
#' Makes a histogram of your data frame, showing the percentage contribution from
#' the top N clones from a chosen dataset in colors.
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param top_clones_choice A data frame of one column, with rownames being the barcodes. Used to get the top clones.
#'@param n_clones Numeric. Number of top clones from each month that should be displayed.
#'@param text_size Numeric. Size of text in plot.
#'@param y_limit Numeric. Limit of y axis in plot.
#'@param other_color Color to use for representation of non-top clones in the bar chart.
#'@param your_title Your title for the plot.
#'@return Displays a bar chart (made by ggplot2) of the samples' top clones in other samples.
#'@examples
#'topclones_barchart(your_data = zh33_Tlineage, top_clones_choice = zh33T1m, y_limit = 50)
#'@export

topclones_barchart <- function(your_data, top_clones_choice, n_clones = 10, text_size = 15,
                               y_limit = 100, other_color = "black", your_title = "Bar Chart"){
  top_barcodes <- rownames(top_clones_choice[order(-top_clones_choice),,drop = FALSE])[1:n_clones]
  your_data <- as.data.frame(100*prop.table(as.matrix(your_data), margin = 2))
  your_data <- rbind(your_data[rownames(your_data) %in% top_barcodes,],
                     colSums(your_data[!(rownames(your_data) %in% top_barcodes),,drop = FALSE]))
  rownames(your_data)[nrow(your_data)] <- "OTHERS"
  your_data <- reshape2::melt(as.matrix(your_data))
  colnames(your_data) <- c("BARCODE", "SAMPLE", "PROP")
  ggplot2::ggplot(your_data, ggplot2::aes(x=SAMPLE, y = PROP,fill = BARCODE))+
    ggplot2::geom_bar(stat="identity", show.legend = FALSE)+
    ggplot2::scale_y_continuous(breaks = seq(0,100,by = 10), name = "Hematopoietic Contribution %")+
    ggplot2::xlab("Sample")+
    ggplot2::theme_classic()+
    ggplot2::ggtitle(your_title)+
    ggplot2::scale_fill_manual(guide = FALSE,
                               values = c(rainbow(length(unique(your_data$BARCODE))-1, s = 1, v = 1, start = 0, end = 0.75, alpha = 1), other_color))+
    ggplot2::theme(text = ggplot2::element_text(size=text_size))+
    ggplot2::coord_cartesian(ylim = c(0, y_limit))


}
