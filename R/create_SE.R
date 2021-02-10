#'@title create_SE
#'
#'@description Creates a SummarizedExperiment object from a data frame containing clonal tracking counts (`your_data`) with rows as observations and columns as samples, and the associated metadata (`meta_data`) with rows as samples and columns of information describing those samples.
#'
#'@param your_data A data frame. For clonal tracking data, this will be individual barcodes or lineage tracing elements in rows and samples in columns.
#'@param meta_data A data frame containing all meta-data. Must, at the very least, include a column called "SAMPLENAME" that contains all of the colnames within the data frame passed as  `your_data` and only those colnames.
#'@param threshold Numeric. The minimum threshold abundance for a barcode to be maintained in the SE. If `threshold_type` is relative, this parameter should be between 0 and 1. If `threshold_type` is absolute, this parameter should be greater than 1. 
#'@param threshold_type Character. One of "relative" or "absolute" relative. If a relative threshold is specified, only those rows which have higher than `threshold` proportion of reads within at least one sample will be kept as non-zero. If an absolute threshold is specified, only those rows which have an absolute read count higher than `threshold` in at least one sample will be kept as non-zero. 
#'@param log_base A numeric indicating which base to use when logging the normalized data
#'@param scale_factor A numeric indicating what scaling factor to use in normalization. For the default value of 1 million, barcode proportions on a per sample basis will be multiplied by 1 million before log+1 normalization.
#'
#'@return Returns a SummarizedExperiment holding your clonal tracking data and the associated metadata.
#'
#'@import SummarizedExperiment
#'
#'@examples
#'create_SE(your_data = read.delim(system.file("sample_data/app_sample_data/sample_data_ZJ31.txt", package = "barcodetrackR"), row.names = 1), meta_data = read.delim(system.file("sample_data/app_sample_data/sample_metadata_ZJ31.txt", package = "barcodetrackR")))
#'
#'@export
#'
create_SE <- function(your_data = NULL,
                      meta_data = NULL,
                      threshold = 0,
                      threshold_type = "relative",
                      log_base = exp(1),
                      scale_factor = 1e6){
  
  # Error checking
  if(is.null(your_data) | is.null(meta_data)){
    stop("NEITHER `your_data` NOR `meta_data` can be NULL")
  }
  if(!setequal(colnames(your_data), meta_data$SAMPLENAME) | ncol(your_data) != nrow(meta_data)){
    stop("The SAMPLENAME column in meta_data must only contain ALL columns in your_data")
  }
  
  if (threshold_type %in% c("relative","absolute") == FALSE){
    stop("The parameter `threshold_type` must be set to `relative` or `absolute`.")
  }
  
  if (threshold_type == "relative"){
    if (threshold < 0 | threshold >= 1){
      stop("Since `threshold_type` is set to `relative`, `threshold` must be greater than or equal to 0 and less than 1.")
    }
  }
  
  if (threshold_type == "absolute"){
    if (threshold <= 1){
      stop("Since `threshold_type` is set to `absolute`, `threshold` must be greater than 1.")
    }
  }
  
  # Thresholding
  pre_thresh_bc_num <- nrow(your_data[rowSums(your_data)>0,])
  if (threshold == 0){
    cat("No threshold supplied. All barcodes will be retained. Be aware that lower abundance barcodes are likely to be less reliable due to sampling bias. To estimate an appropriate threshold, please see the barcodetrackR function `estimate_barcode_threshold`. \n")
  }
  else if(threshold != 0){
    your_data <- threshold(your_data, thresh = threshold, thresh_type = threshold_type)
    cat("Removed", pre_thresh_bc_num - nrow(your_data), "barcodes from the supplied dataframe based on", threshold_type, "threshold of", threshold, "\n")
  }
  
  # Check that all samples contain data
  if(any(colSums(your_data) == 0)){
    bad_samples <- which(colSums(your_data) == 0)
    if (threshold == 0){
      cat("The following samples have no data. \n")
      cat(colnames(your_data)[bad_samples], sep = ", ")
      cat("\n")
      stop("Please remove the empty sample before instantiating the SE.")
    } else if (threshold != 0){
      cat("The following samples have no data after thresholding. \n")
      cat(colnames(your_data)[bad_samples], sep = ", ")
      at("\n")
      stop("Please try a more permissive threshold or remove the sample(s) prior to instantiating the SE.")
    }
  }
  
  your_data.ranks <-  as.data.frame(apply(-your_data, 2, rank, ties.method = "min", na.last = "keep"))
  your_data.proportions <-  as.data.frame(prop.table(as.matrix(your_data),2))
  your_data.normalized <- your_data.proportions * scale_factor
  your_data.logged <- log(1+your_data.normalized, base = log_base)
  your_SE <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = your_data,
                                                                      proportions = your_data.proportions,
                                                                      ranks = your_data.ranks,
                                                                      normalized = your_data.normalized,
                                                                      logs = your_data.logged),
                                                        colData=meta_data)

  S4Vectors::metadata(your_SE)$scale_factor <- scale_factor
  S4Vectors::metadata(your_SE)$log_base <- log_base
  if (threshold == 0){
    S4Vectors::metadata(your_SE)$threshold_type <- "none"
    S4Vectors::metadata(your_SE)$threshold_value <- threshold
  } else if (threshold != 0){
    S4Vectors::metadata(your_SE)$threshold_type <- threshold_type
    S4Vectors::metadata(your_SE)$threshold_value <- threshold
  }
  
  return(your_SE)
}

