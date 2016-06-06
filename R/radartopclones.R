#' Radar Top Clones
#'
#' Makes a radar chart of your data frame, plotting the percentage contribution from
#' the top N clones from a specified column.
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param columnChoice_Name Character. The name of the column to take the n clones from.
#'@param n_clones Numeric. The number of top clones to take from columnChoice_Name.
#'@return labelsize Numeric. Size of the labels.
#'@examples
#'
#'radartopclones(your_data = zh33[,c(1:6)], columnChoice_Name = "zh331mT.fastq",
#'               your_title = "1m Radar Plot")
#'@export


radartopclones <- function(your_data, columnChoice_Name = "", n_clones = 10, your_title = "", labelsize = 1){
  maxRow = rep(1,ncol(your_data))
  minRow = rep(0,ncol(your_data))
  your_data <- prop.table(as.matrix(your_data), margin = 2)
  #takes care of limiting top clones to the number of barcodes in a sample
  if (n_clones > min(colSums(your_data!=0)))
    stop("Number of clones exceeds number of barcodes for a lineage")
  dfbarcodes <- sum(your_data[,columnChoice_Name]!=0)
  if (n_clones > dfbarcodes)
    n_clones <- dfbarcodes
  your_data <- your_data[order(-your_data[,columnChoice_Name]),][c(1:n_clones),]
  if (n_clones > 1)
    your_data <- colSums(your_data)
  your_data <- rbind(maxRow, minRow, your_data)
  fmsb::radarchart(as.data.frame(your_data), title = your_title, centerzero = TRUE, vlcex = labelsize)
}
