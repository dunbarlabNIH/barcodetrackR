#'@title Binary Barcode Heat Map
#'
#'@description Creates a binary heatmap showing the presence of new clones from L to R in the dataset.
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param your_names Names for each column of the data frame.
#'@param threshold A chosen threshold that decides if a barcode exists or not.
#'@param grid Logical. Include a grid or not in the heatmap.
#'@param col_labels The size of the column labels.
#'@return Displays a binary heat map in the current plot window.
#'@examples
#'BBHM(your_data = zh33, names = colnames(zh33), threshold = 1000,
#'       grid = TRUE, col_labels = 3)
#'BBHM(your_data = zh33, grid = FALSE)
#'@export


BBHM <- function(your_data, your_names = colnames(your_data), threshold = 0,
                 grid = TRUE, col_labels = 1){
  your_data <- your_data[,your_names]
  your_data <- your_data[rowSums(your_data) != 0,]
  your_data[your_data <= threshold] = 0
  your_data[your_data > threshold] = 1
  your_data <- your_data[do.call(order, -as.data.frame(your_data)),]
  gplots::heatmap.2(as.matrix(your_data),
            scale = "none",
            dendrogram = "none",
            Colv = FALSE,
            labRow = "",
            col = colorRampPalette(c("white", "blue"))(2),
            key = FALSE,
            density.info = "none",
            trace = "none",
            Rowv = FALSE,
            colsep = if (grid) 1:ncol(your_data) else NULL,
            rowsep = if (grid) 1:nrow(your_data) else NULL,
            cexCol = col_labels)
}
