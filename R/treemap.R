#' Barcode treemap
#'
#' Makes a treemap of the chosen column in your data frame.
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param column_choice Character. The name of the column to make the treemap from.
#'@param threshold Numeric. Threshold for the whole data frame. Anything lower than this becomes 0.
#'@param colors_choice Character vector of size 2 specifying which colors to use for the treemap range.
#'@param your_title Character. Title for the plot.
#'@return Displays a treemap (made by treemap) of the chosen column.
#'@examples
#'barcode_treemap(your_data = zh33, column_choice = "32m_CD4.fastq", colors_choice = c("red", "green"))
#'@export
barcode_treemap <- function(your_data, column_choice = NULL, threshold = 0, colors_choice = c("blue", "yellow"), your_title = ""){
  your_data[your_data < threshold] <- 0
  your_data <- as.data.frame(prop.table(as.matrix(your_data), margin = 2))
  your_data <- your_data[,column_choice, drop = FALSE]
  your_data <- as.data.frame(melt(as.matrix(your_data)))
  treemap::treemap(your_data, index = "Var1", vSize = "value", type = "manual", title = your_title, palette = colorRampPalette(colors_choice)(10), fontsize.labels = 0, fontsize.title = 30)
}
