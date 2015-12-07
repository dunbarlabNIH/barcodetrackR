#'@title Clonal Diversity
#'
#'@description Gives the diversity for each sample (column) in your data.
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param threshold Numeric. Universal threshold applied to your_data (anything below this threshold will be set to 0.#'@param writekey Logical. Write contributions to file or not.
#'@param divindex Character. One of "shannon", "simpson", or "invsimpson".
#'@return Vector of the diversities of each column.
#'@examples
#'clonaldiversity(your_data = zh33, threshold = 2000, divindex = "shannon")
#'

clonaldiversity <- function(your_data, threshold = 2000, divindex = "shannon") {
  your_data[your_data < threshold] <- 0
  return(vegan::diversity(your_data, index = divindex, MARGIN = 2))
}
