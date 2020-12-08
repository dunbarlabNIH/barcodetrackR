#'@title create_SE
#'
#'@description Creates a SummarizedExperiment object from a data frame containing clonal tracking counts (`your_data`) with rows as observations and columns as samples, and the associated metadata (`meta_data`) with rows as samples and columns of information describing those samples.
#'
#'@param your_data A data frame. For clonal tracking data, this will be individual barcodes or lineage tracing elements in rows and samples in columns.
#'@param meta_data A data frame containing all meta-data. Must, at the very least, include a column called "SAMPLENAME" that contains all of the colnames within the data frame passed as  `your_data` and only those colnames.
#'@param threshold Optional. A numeric indicating which reads to filter out (only those rows which have higher than `threshold` proportion of reads within at least one sample will be kept as non-zero) prior to creating the SE object.
#'@param log_base A numeric indicating which base to use when logging the percentages.
#'@param scale_factor A numeric indicating what scaling factor to use in normalization. For the default value of 1 million, barcode percentages on a per sample basis will be multiplied by 1 million before log+1 normalization.
#'
#'@return Returns a SummarizedExperiment holding your clonal tracking data and the associated metadata.
#'
#'@import SummarizedExperiment
#'
#'@examples
#'create_SE(your_data = exp1_clonal_counts_data, meta_data = exp1_meta_data)
#'
#'@export
#'
create_SE <- function(your_data = NULL,
                      meta_data = NULL,
                      threshold = 0,
                      log_base = exp(1),
                      scale_factor = 1e6){
  if(is.null(your_data) | is.null(meta_data)){
    stop("NEITHER `your_data` NOR `meta_data` can be NULL")
  }
  if(!setequal(colnames(your_data), meta_data$SAMPLENAME) | ncol(your_data) != nrow(meta_data)){
    stop("The SAMPLENAME column in meta_data must only contain ALL columns in your_data")
  }
  if(threshold > 0){
    your_data <- threshold(your_data, thresh = threshold)
  }
  if(any(colSums(your_data) == 0)){
    stop("One of your samples has no data after thresholding")
  }
  your_data.ranks <-  as.data.frame(apply(-your_data, 2, rank, ties.method = "min", na.last = "keep"))
  your_data.percentages <-  as.data.frame(prop.table(as.matrix(your_data),2))
  your_data.normalized <- your_data.percentages * scale_factor
  your_data.logged <- log(1+your_data.normalized, base = log_base)
  your_SE <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = your_data,
                                                                      percentages = your_data.percentages,
                                                                      ranks = your_data.ranks,
                                                                      normalized = your_data.normalized,
                                                                      logs = your_data.logged),
                                                        colData=meta_data)

  S4Vectors::metadata(your_SE)$scale_factor <- scale_factor
  S4Vectors::metadata(your_SE)$log_base <- log_base


  return(your_SE)
}

