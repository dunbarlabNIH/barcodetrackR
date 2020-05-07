#' Stat histogram
#'
#' Given a summarized experiment, gives a histogram of any barcoding assay or choice of metadata. 
# '
#'@param your_SE Your SummarizedExperiment of barcode data and associated metadata
#'@param data_choice Either "barcode_stats" which allows you to see the distribution of any assay for chosen samples or "aggregate_stats" which allows you to see distirbution of a piece of metadata cross all samples.
#'@param assay_choice When data_choice is set to "barcode_stats", provide which assay fo the SE you want to display in the histogram.
#'@param sample_select When data_choice is set to "barcode_stats", provide a list of sample names to display histograms for
#'@param metadata_stat When data_choice is set to "aggregate_stats", The column of metadata to display in a histogram. Must be numeric.
#'@param group_by When data_choice is set to "aggregate_stats", Group and color the histogram by another column of metadate. If NULL, no grouping applied
#'@param scale_all_y Logical. Whether or not to plot all plots on the same y axis limits.
#'@param y_log_axis Logical. Whether or not to put y axis on log scale
#'@param x_log_axis Logical. Whether or not to put x axis on log scale
#'@param n_bins Number of bins for histogram. Default is 10.
#'@param n_cols Number of columns for faceted histograms.
#'@param text_size Size of text.
#'@param alpha The transparency of the histograms if group_by is provided. Lower = more transparent. 1 = opaque
#'@return Histogram of chosen statistics
#'
#'@importFrom rlang %||%
#'@examples
#'stat_hist(your_se = SE, metadata_stat = "GFP_percent")
#'@export
stat_hist <- function(your_SE,
                      data_choice = "barcode_stats",
                      assay_choice = "counts",
                      sample_select = NULL,
                      metadata_stat = NULL,
                      group_by = NULL,
                      scale_all_y = FALSE,
                      y_log_axis = FALSE,
                      text_size = 12,
                      n_bins = 10,
                      n_cols = 1,
                      alpha = 0.5){
  
  
  
  if (data_choice == "barcode_stats"){
    
    # Error handling
    if (assay_choice %in% names(assays(your_SE)) == FALSE){
      stop("chosen assay is not in your_SE.")
    }
    
    # Load data
    your_data <- SummarizedExperiment::assays(your_SE)[[assay_choice]]
    
    # sample selection
    sample_select <- sample_select %||% colnames(your_SE)[1]
    
    # Error handling
    if (any(sample_select %in% colnames(your_SE) == FALSE)){
      stop("one of the chosen samples is not a sample name in your_SE.")
    }
    
    # Subset the data
    your_data = as.data.frame(your_data[,sample_select])
    
    # Make sure the order is correct
    
  } else if (data_choice == "aggregate_stats"){
    
    # Load metadata
    meta_data <- as.data.frame(SummarizedExperiment::colData(your_SE))
    
    # Error handling
    if (metadata_stat %in% colnames(meta_data) == FALSE){
      stop("metadata_stat is not a piece of colData.")
    }
    
    if(class(meta_data[,metadata_stat]) != "numeric"){
      stop("metadata_stat is not a numeric vector.")
    } 
    
    if (is.null(group_by) == FALSE){
      if (group_by %in% colnames(meta_data) == FALSE){
        stop("group_by is not a piece of colData.")
      }
    }
    
  } else {
    stop("data_choice must be one of 'barcode_stats' or 'aggregate_stats' .")
  }
  
  # Choice for y and x axis
  y_axis_choice <-"identity"
  if(y_log_axis){
    y_axis_choice <- "pseudo_log"
  }

  
  # Plotting barcode_stats
  if (data_choice == "barcode_stats"){
    
    plot_list <- lapply(1:length(sample_select), function(i){
    
      g <- suppressMessages(
        ggplot2::ggplot(your_data, ggplot2::aes(x = your_data[,i]))+
        ggplot2::geom_histogram(bins = n_bins, color = "white", fill = "dodgerblue2")+
        ggplot2::theme_classic()+
        ggplot2::labs(x = paste0("barcode ",assay_choice), title = sample_select[i])+
        ggplot2::theme(text = ggplot2::element_text(size = text_size))+
        ggplot2::scale_y_continuous(trans = y_axis_choice)
      )
      
    })
    
    if(scale_all_y){
      plot_max <- max(unlist(lapply(plot_list, function(x){ggplot2::ggplot_build(x)$layout$panel_params[[1]]$y.range[2]})))
      plot_list <- lapply(plot_list, function(x){x + ggplot2::coord_cartesian(ylim = c(0, plot_max))})
    }
    
    p <- cowplot::plot_grid(plotlist = plot_list, ncol = n_cols)
  
    }
    
  
  # Plotting aggregate_stats
  if (data_choice == "aggregate_stats"){
  
    # Grouped histogram
    if (is.null(group_by) == FALSE){
      p <- suppressWarnings(
        ggplot2::ggplot(meta_data, ggplot2::aes(x = meta_data[,metadata_stat], color = meta_data[,group_by], fill = meta_data[,group_by]))+
        ggplot2::geom_histogram(bins = n_bins, alpha = alpha, position = "identity")+
        ggplot2::theme_classic()+
        ggplot2::labs(x = paste0("metadata: ",metadata_stat), fill = group_by, title = "Aggregate Statistics") +
        ggplot2::guides(color = FALSE)+
        ggplot2::theme(text = ggplot2::element_text(size = text_size))+
        ggplot2::scale_y_continuous(trans = y_axis_choice)
      )
      
    } else {   # Regular histogram
      p <- suppressWarnings(
        ggplot2::ggplot(meta_data, ggplot2::aes(x = meta_data[,metadata_stat]))+
          ggplot2::geom_histogram(bins = n_bins, color = "white", fill = "dodgerblue2")+
          ggplot2::theme_classic()+
          ggplot2::labs(x = paste0("metadata: ",metadata_stat), title = "Aggregate Statistics")+
          ggplot2::theme(text = ggplot2::element_text(size = text_size))+
          ggplot2::scale_y_continuous(trans = y_axis_choice)
      )
      
    }
    
    
  }

  p
  
}
