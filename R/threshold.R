#' Threshold
#'
#' This function takes in sequence data in table form, along with a threshold,
#" and returns new table with the criteria that at least one cell type is above
#' the threshold.
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param thresh Numeric.
#'@return A data frame where all rows (barcodes) that did not have at least one elemnt meet the threshold have been discarded.
#'@examples
#'threshold(zh33, thresh = 2000)
#'@export



threshold = function(your_data, thresh=2000) {

#makes vector of thresholds for each fq file, assuming 4 million reads in original fq file
#note, threshold given in arguments is scaled to number of actual reads...
#here it's 2000 reads that should be present in 4000000, which is then made into a percentage and
#scaled to the number of reads in a file
#note that 2,000/4,000,000 is 0.0005, or 0.05%
threshes <- (colSums(your_data)*(thresh/4000000))

#keeps rows that have at least one element larger than thresh for that column
thresholded_data <- your_data[apply(your_data, 1, function(x) {any(threshes < x)}),]

#in the unlikely case that a file has 0 reads, this will prevent the phantom
#column from appearing in the final data frame
thresholded_data[is.na(thresholded_data)] <- 0

#returns subsetted data
return(thresholded_data);

}
