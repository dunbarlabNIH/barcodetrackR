#'@title Top Clones Getter
#'
#'@description Given a subset of data, returns only the rows that contain the top N entries of each column.
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param top Numeric. How many top entries to find.
#'@return Vector of indices of the top N entries of each column.
#'@examples
#'gettopindices(your_data = zh33[,c(1:5)], top = 10)
#'@export

gettopindices <- function(your_data, top = 5){
  order.list <- apply(your_data, 2, function(x){sort.list(-x)[1:top]})
  top_clones <- unique(as.vector(order.list))
  return(top_clones)
}
