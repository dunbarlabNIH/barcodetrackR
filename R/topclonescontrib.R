#' Top Clones Contributions
#'
#' Makes a line plot of your data frame, plotting the percentage contribution from
#' the top N clones from each column as a line.
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param n_clones Numeric. The number of top clones to take from columnChoice_Name.
#'@param linesize Numeric. Size of the lines.
#'@param pointsize Numeric. Size of the points.
#'@param your_title Character. Title for the plot.
#'@return Prints a plot (using ggplot2) of the percentage contribution from each sample.
#'@examples
#'topclonescontrib(your_data = zh33[,c(1:6)], n_clones = 100, your_title = "SAMPLE")
#'@export

topclonescontrib <- function(your_data, n_clones = 10, linesize = 2, pointsize = 3, your_title = ""){

  your_data <- as.data.frame(100*prop.table(as.matrix(your_data),2))

  #takes care of limiting top clones to the number of barcodes in a sample
  if (n_clones > min(colSums(your_data!=0)))
    stop("Number of clones exceeds number of barcodes for a lineage")

  newdf <- matrix(ncol = ncol(your_data), nrow = ncol(your_data))
  colnames(newdf) <- colnames(your_data)

  for(i in 1:ncol(your_data)){
    temp <- your_data[order(-your_data[,i]),][c(1:n_clones),]
    if (n_clones > 1)
      temp <- colSums(temp)
    newdf[i,] <- temp
  }

  rownames(newdf) <- paste0("top_",n_clones,"_",colnames(your_data))

  newdf <- reshape2::melt(newdf)

  ggplot2::ggplot(newdf, ggplot2::aes(x=Var2, y= value, colour = Var1, group = Var1))+
    ggplot2::geom_line(size = linesize)+
    ggplot2::geom_point(fill = "white", size = pointsize)+
    ggplot2::coord_cartesian(ylim = c(0,1))+
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.1))+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1))+
    ggplot2::ylab("Hematopoietic Contribution")+
    ggplot2::xlab("")+
    ggplot2::ggtitle(your_title)

}
