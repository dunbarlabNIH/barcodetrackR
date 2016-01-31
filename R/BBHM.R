#'@title Binary Barcode Heat Map
#'
#'@description Creates a binary heatmap showing the presence of new clones from L to R in the dataset.
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param threshold A chosen threshold that decides if a barcode exists or not.
#'@param col_labels The size of the column labels.
#'@return Displays a binary heat map in the current plot window.
#'@examples
#'BBHM(your_data = zh33, threshold = 1000, col_labels = 3)
#'BBHM(your_data = zh33)
#'@export


BBHM <- function(your_data, threshold = 0, col_labels = 1){
  your_data[your_data <= threshold] = 0
  your_data[your_data > threshold] = 1
  your_data <- your_data[rowSums(your_data) != 0,]
  if(any(colSums(your_data) == 0)){
    stop("BBHM: One of your columns had no reads above the threshold.")
  }
  your_data <- your_data[do.call(order, -as.data.frame(your_data)),]
  print(head(your_data))
  gplots::heatmap.2(as.matrix(your_data),
            scale = "none",
            dendrogram = "none",
            Colv = FALSE,
            labRow = "",
            col = colorRampPalette(c("white", "navyblue"))(2),
            key = FALSE,
            density.info = "none",
            trace = "none",
            Rowv = FALSE,
            sepcolor = "black",
            colsep = c(0, ncol(your_data)),
            rowsep = c(0, nrow(your_data)),
            sepwidth = c(1/ncol(your_data)^4, 1/nrow(your_data)^3),
            srtCol = 45,
            lhei = c(1,10),
            lwid = c(1,11),
            margins = c(10,9),
            cexCol = col_labels
            )

}
