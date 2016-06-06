#' Data Frame Merger
#'
#' Merges an arbitrary number of data frames together by row name, including all entries. Hi, mom!
#'@param ... A list of data frames to merge by row name.
#'@return Merged data frame.
#'@examples
#'merger(dataframe1, dataframe2, dataframe3)
#'@export


merger <- function(...){
  for(i in list(...)){
    if(!is.data.frame(i) && !is.vector(i))
      stop("All arguments must be data frames or vectors!")
  }
  your_data <- data.frame()
  for(i in list(...)){
    if(is.vector(i)){
      i <- as.data.frame(i)
    }
    your_data <- merge(your_data,i, by = 0, all = T)
    rownames(your_data) <- your_data$Row.names
    your_data$Row.names <- NULL
  }
  return(your_data)
}
