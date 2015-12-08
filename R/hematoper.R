#' Hematopoietic Contribution
#'
#' Usually used for tracking a cell lineage's top clones over time. Best to
#' put samples in order from left to right (i.e. "1mT", "2mT", "3mT", etc.)
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param months Numeric vector. Vector corresponding to the months of each sample.
#'@param n_clones Numeric. Number of top clones to view from each sample.
#'@param scale_percent Logical. Scale as percent of hematopoiesis or keep as number of reads.
#'@param linesize Numeric. Thickness of the lines.
#'@return Displays a stacked area plot (made by ggplot2) of the samples' top clones.
#'@examples
#'hematoper(your_data = zh33_Tlineage, months = c(1,2,3,4,5,6.5, 9.5, 12), scale_percent = TRUE)
#'@export

hematoper <- function(your_data, months, n_clones = 10, scale_percent = FALSE, linesize = 0.5){
  if(scale_percent == TRUE) {
    your_data <- prop.table(as.matrix(your_data), margin = 2)
  }
  your_data <- your_data[barcodetrackR::gettopindices(your_data, top = n_clones),]
  colnames(your_data) <- months
  melty <- reshape2::melt(as.matrix(your_data))
  ggplot2::ggplot(melty, ggplot2::aes(x=Var2, y = value, group = Var1, fill = Var1))+
    ggplot2::geom_area(position = "stack", linetype = 1, size = linesize, colour = "black")+
    ggplot2::theme(legend.position = "none")+
    ggplot2::ylab("Hematopoietic Contribution")+
    ggplot2::xlab("")
}
