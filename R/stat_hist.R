#' Stat histogram
#'
#' Given a summarized experiment, gives a histogram of any piece of metadata. 
# '
#'@param your_SE Your SummarizedExperiment of barcode data and associated metadata
#'@param metadata_stat The column of metadata to display in a histogram. Must be numeric.
#'@param group_by Group and color the histogram by another column of metadat. If NULL, no grouping applied
#'@param text_size Size of text.
#'@param alpha The transparency of the histograms if group_by is provided. Lower = more transparent. 1 = opaque
#'@examples
#'stat_hist(your_se = SE, metadata_stat = "GFP_percent")
#'@export
stat_hist <- function(your_SE, metadata_stat, group_by = NULL, text_size = 20, alpha = 0.5){
  
  # Load data
  # your_data <- SummarizedExperiment::assays(your_SE)$counts
  meta_data <- as.data.frame(SummarizedExperiment::colData(your_SE))
  
  # Basic error handling
  if (metadata_stat %in% colnames(meta_data) == FALSE){
    stop("metadata_stat is not a piece of colData.")
  }
  
  if (is.null(group_by) == FALSE){
    if (group_by %in% colnames(meta_data) == FALSE){
      stop("group_by is not a piece of colData.")
    }
  }

  if(class(meta_data[,metadata_stat]) == "numeric"){
    stop("metadata_stat is not a numeric vector.")
  }
  
  # Grouped histogram
  if (is.null(group_by) == FALSE){
    p <- ggplot2::ggplot(meta_data, aes(x = meta_data[,metadata_stat], color = meta_data[,group_by], fill = meta_data[,group_by]))+
      ggplot2::geom_histogram(alpha = alpha, position = "identity")+
      ggplot2::theme_classic()+
      ggplot2::labs(x = metadata_stat, fill = group_by) +
      ggplot2::guides(color = FALSE)+
      ggplot2::theme(text = ggplot2::element_text(size = text_size))
    
  } else {   # Regular histogram
    p <- ggplot2::ggplot(meta_data, aes(x = meta_data[,metadata_stat]))+
      ggplot2::geom_histogram(color = "darkblue", fill = "lightblue")+
      ggplot2::theme_classic()+
      ggplot2::labs(x = metadata_stat)+
      # ggplot2::guides(fill= FALSE, color = FALSE)+
      ggplot2::theme(text = ggplot2::element_text(size = text_size))
    
  }

  p
  
}
