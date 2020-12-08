#'@title subset_SE
#'
#'@description Subsets an existing SummarizedExperiment object.
#'
#'@param your_SE A SummarizedExperiment object.
#'@param ... Arguments passed to subset_SE in the form of `X = keys` where `X` is a column from SE's colData and `keys` are entries in the colData to subset.
#'@return Returns a subsetted SummarizedExperiment object.
#'@examples
#'zh33 <- create_SE(your_data = ZH33_file, meta_data = ZH33_meta_data)
#'zh33.B.2month <- subset_SE(zh33, Cell_type = "B", Timepoint = "2m")
#'@export
subset_SE <- function(your_SE, ...) {
  arguments <- list(...)
  if(length(arguments) > 0){
    for(i in 1:length(arguments)){
      subset_name <- names(arguments)[[i]]
      subset_vars <- arguments[[i]]
      your_SE <- your_SE[, your_SE[[subset_name]] %in% subset_vars]
    }
  }
  return(your_SE)
}
