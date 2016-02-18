#' Rank Abundance
#'
#' Rank-Abundance Plot of the Barcodes for each sample chosen. Uses the cumulative sum so that all
#' samples eventually reach 100%.
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param dot_size Numeric. Size of the points.
#'@param text_size Numeric. Size of text in plot.
#'@return Displays a rank-abundance plot (made by ggplot2) of the samples chosen.
#'@examples
#'rank_abundance_plot(your_data = zh33_Tlineage[,c(1,2,4)], dot_size = 4)
#'@export
rank_abundance_plot <- function(your_data, dot_size = 5, text_size = 10){
  your_data <- your_data[rowSums(your_data) != 0,]
  your_data <- 100*prop.table(as.matrix(your_data), margin = 2)
  for(i in 1:ncol(your_data)){
    your_data[,i] <- your_data[order(-your_data[,i]),i]
    your_data[,i] <- cumsum(your_data[,i])
  }
  your_data <- reshape2::melt(your_data)
  colnames(your_data) <- c("Barcode", "Sample", "Prop")
  your_data$Rank <- 0
  samples <- unique(your_data$Sample)
  for(x in samples){
    your_data[your_data$Sample == x,]$Rank <- rank(your_data[your_data$Sample == x,]$Prop, ties.method = "min")
  }
  ggplot2::ggplot(data = your_data, ggplot2::aes(x = Rank, y = Prop, group = Sample))+
    ggplot2::geom_point(ggplot2::aes(fill = Sample, colour = Sample), shape = 21, size = dot_size)+
    ggplot2::theme_bw()+
    ggplot2::scale_x_continuous(breaks = seq(0,max(your_data$Rank),by = 200))+
    ggplot2::scale_y_continuous(breaks = seq(0,100,by = 10))+
    ggplot2::theme(text = ggplot2::element_text(size = text_size))
}
