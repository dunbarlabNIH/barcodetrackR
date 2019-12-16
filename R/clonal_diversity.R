#'@title Clonal Diversity
#'
#'@description Gives the diversity for each sample (column) in your data.
#'
#'@param your_data A data frame. Individual barcodes in rows and samples in columns.
#'@param divindex Character. One of "shannon", "simpson", or "invsimpson".
#'@return Vector of the diversities of each column.
#'@examples
#'clonaldiversity(your_data = zh33, divindex = "shannon")
#'@export
clonal_diversity <- function(your_data, threshold = 2000, divindex = "shannon") {
  your_data[your_data < threshold] <- 0
  return(vegan::diversity(your_data, index = divindex, MARGIN = 2))
}
