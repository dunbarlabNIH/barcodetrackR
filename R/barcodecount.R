#' Barcode Counter
#'
#' Counts the number of barcodes sequentially by column cumulatively, uniquely, or by first presence (new).
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param count Character. One of "cumulative", "unique", or "new".
#'@param months Numeric vector. Corresponds to months of each column of your_data.
#'@param threshold Numeric. Universal threshold applied to your_data (anything below this threshold will be set to 0.
#'@return Vector containing the barcode count of each column, named by month.
#'@examples
#'barcodecount(your_data = zh33, count = "unique", months = c(1,2,3,4.5,6.5,9), threshold = 1000)
#'barcodecount(your_data = zg66, count = "new", months = c(1:10), threshold = 900)
#'
#'
#'
barcodecount <- function(your_data, count = "cumulative", months = NULL, threshold = 0) {
  if (length(months) != ncol(your_data))
    stop("Number of months should match number of columns")
  if (!(count %in% c("cumulative", "unique", "new")))
    stop("Count must be one of \"cumulative\", \"unique\", or \"new\" ")
  your_data[your_data < threshold] <- 0
  if (count == "unique"){
    vec <- colSums(your_data > 0)
  }
  if (count == "cumulative"){
    vec <- colSums(t(apply(your_data > 0, 1, cumsum)) > 0)
  }
  if (count == "new"){
    vec <- colSums(t(apply(t(apply(your_data > 0, 1, cumsum)) > 0,1, cumsum)) == 1 )
  }
  names(vec) <- months
  return(vec)
}
