#' Threshold
#'
#' This function takes in sequence data in table form, along with a threshold,
#" and returns new table with the criteria that each row contains at least one element
#" that registers above the threshold. Threshold here is defined as fractional contribution
#' to each column (e.g. if threshold is set as 0.0005, only rows in which an element is above 0.05% of
#' its column will be kept).
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param thresh Numeric.
#'@return A data frame where all rows (barcodes) that did not have at least one element meet the threshold have been discarded.
#'@examples
#'threshold(zh33, thresh = 0.0005)
#'@export

threshold = function(your_data, thresh=0.0005) {

  #makes vector of thresholds for each fq file, assuming 4 million reads in original fq file
  #note that s 0.0005, or 0.05%
  threshes <- (colSums(your_data)*(thresh))

  #keeps rows that have at least one element larger than thresh for that column
  thresholded_data <- your_data[apply(your_data, 1, function(x) {any(threshes < x)}),]

  #in the unlikely case that a file has 0 reads, this will prevent the phantom
  #column from appearing in the final data frame
  thresholded_data[is.na(thresholded_data)] <- 0

  #returns subsetted data
  return(thresholded_data)

}
