#'@title Threshold SE
#'
#'@description Removes barcodes from a SummarizedExperiment object which have an abundance lower than the provided relative or absolute threshold. See the function `estimate_barcode_threshold` to estimate an appropriate threshold for an SE.
#'
#'@param your_SE A Summarized Experiment object.
#'@param threshold_value Numeric. The minimum threshold abundance for a barcode to be maintained in the SE. If `threshold_type` is relative, this parameter should be between 0 and 1. If `threshold_type` is absolute, this parameter should be greater than 1. 
#'@param threshold_type Character. One of "relative" or "absolute" relative. If a relative threshold is specified, only those rows which have higher than `threshold_value` proportion of reads within at least one sample will be kept as non-zero. If an absolute threshold is specified, only those rows which have an absolute read count higher than `threshold_value` in at least one sample will be kept as non-zero. 
#'@param verbose Logical. If TRUE, print the total number of barcodes removed from the SE.
#'
#'@return Returns a SummarizedExperiment containing only barcodes which passed the supplied threshold in at least one sample. All of the defualt assays are re-calculated after thresholding is applied. Note that since tthe SE is re-instantiated, any custom assays should be recalculated after thresholding.
#'
#'@import SummarizedExperiment
#'
#'@examples
#'threshold_SE(your_SE = wu_subset, threshold_value = 0.005, threshold_type = "relative", verbose = TRUE)
#'
#'@export
#'
threshold_SE <- function(your_SE,
                         threshold_value,
                         threshold_type = "relative",
                         verbose = TRUE){
  
  # Error checking
  if (threshold_type %in% c("relative","absolute") == FALSE){
    stop("The parameter `threshold_type` must be set to `relative` or `absolute`.")
  }
  
  if (threshold_type == "relative"){
    if (threshold_value <= 0 | threshold_value >= 1){
      stop("Since `threshold_type` is set to `relative`, `threshold_value` must be greater than 0 and less than 1.")
    }
  }
  
  if (threshold_type == "absolute"){
    if (threshold_value <= 1){
      stop("Since `threshold_type` is set to `absolute`, `threshold_value` must be greater than 1.")
    }
  }
  
  # Get number of barcodes before thresholding
  pre_thresh_bc_num <- nrow(SummarizedExperiment::assays(your_SE)$counts)
  
  # Apply threshold
  your_data <- threshold(SummarizedExperiment::assays(your_SE)$counts, thresh = threshold_value, thresh_type = threshold_type)
  
  # Check that all samples still contain data
  if(any(colSums(your_data) == 0)){
    bad_samples <- which(colSums(your_data) == 0)
    cat("The following samples have no data after thresholding. \n")
    cat(colnames(SummarizedExperiment::assays(your_SE)$counts)[bad_samples], sep = ", ")
    cat("\n")
    stop("Please try a more permissive threshold or remove the sample(s) prior to thresholding")
  }
  
  # Re-calculate other assays
  your_data.ranks <-  as.data.frame(apply(-your_data, 2, rank, ties.method = "min", na.last = "keep"))
  your_data.proportions <-  as.data.frame(prop.table(as.matrix(your_data),2))
  your_data.normalized <- your_data.proportions * your_SE@metadata$scale_factor
  your_data.logged <- log(1+your_data.normalized, base = your_SE@metadata$log_base)
  
  # Create new thresholded SE
  thresh_SE <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = your_data,
                                                                        proportions = your_data.proportions,
                                                                        ranks = your_data.ranks,
                                                                        normalized = your_data.normalized,
                                                                        logs = your_data.logged),
                                                          colData=your_SE@colData)
  
  # Update metadata with thresholding information
  S4Vectors::metadata(thresh_SE) <- S4Vectors::metadata(your_SE)
  S4Vectors::metadata(thresh_SE)$threshold_type <- threshold_type
  S4Vectors::metadata(thresh_SE)$threshold_value <- threshold_value
  
  # Print number of barcodes removed
  if (verbose){
    cat("Removed", pre_thresh_bc_num - nrow(SummarizedExperiment::assays(thresh_SE)$counts), "barcodes from the supplied dataframe based on", threshold_type, "threshold of", threshold_value, "\n")
  }
  
  return(thresh_SE)
}

