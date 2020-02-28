#' Ridge plot
#'
#' Given a summarized experiment, gives ridge plots showing percent total contribution to both lineages.
# '
#'@param your_SE Your SummarizedExperiment of barcode data and associated metadata
#'@param split_bias_on The column of metadata corresponding to cell types
#'@param bias_1 The first cell type to be compared. Will be on the RIGHT side of the ridge plot
#'@param bias_2 The second cell type to be compared. Will be on the LEFT side of the ridge plot
#'@param split_bias_over The column of metadata to plot by. If numeric, y axis will be in increasing order. If categorical, it will follow order of metadata.
#'@param bias_over Choice(s) from the column designated in `split_bias_over` that will be used for plotting. Defaults to all.
#'@param weighted If true, the density estimation will be weighted by the overall contribution of each barcode
#'@param text_size The size of the text in the plot.
#'@param scale Sets the overlap between ridges. Larger values give more overlap. 
#'@return Bias plot for two lineages over time.
#'
#'@importFrom dplyr rename
#'@import ggridges
#'@import ggplot2
#'
#'@examples
#'ridge_plot(your_se = SE, split_bias_on = "selection_type", bias_1 = "B", bias_2 = "T", split_bias_over = "Timepoint")
#'@export
ridge_plot <- function(your_SE,
                       split_bias_on,
                       bias_1,
                       bias_2,
                       split_bias_over,
                       bias_over = NULL,
                       weighted = F,
                       text_size = 16,
                       scale = 1){
  
  # Load data
  your_data <- SummarizedExperiment::assays(your_SE)$counts
  meta_data <- SummarizedExperiment::colData(your_SE)

  # Basic error handling
  coldata_names <- colnames(SummarizedExperiment::colData(your_SE))
  if(any(! c(split_bias_on) %in% coldata_names)){
    stop("split_bias_on must match a column name in colData(your_SE)")
  }
  if(! bias_1 %in% unique(SummarizedExperiment::colData(your_SE)[[split_bias_on]])){
    stop("bias_1 is not in split_bias_on")
  }
  if(! bias_2 %in% unique(SummarizedExperiment::colData(your_SE)[[split_bias_on]])){
    stop("bias_2 is not in split_bias_on")
  }
  
  # Only keep data that matches filter
  your_data <- your_data[,meta_data[,split_bias_on] == bias_1 | meta_data[,split_bias_on] == bias_2]
  meta_data <- meta_data[meta_data[,split_bias_on] == bias_1 | meta_data[,split_bias_on] == bias_2,]
  
  # Make the split_bias_on an ordered factor
  meta_data[,split_bias_on] <- factor(meta_data[,split_bias_on], levels = c(bias_1,bias_2))
  
  # error handling
  if(any(! c(split_bias_over) %in% coldata_names)){
    stop("split_bias_over must match a column name in colData(your_SE)")
  }
  
  if(is.numeric(SummarizedExperiment::colData(your_SE)[[split_bias_over]])){
    # Sort
    bias_over <- bias_over %||% sort(meta_data[,split_bias_over])
    # Change to factor
    bias_over <- factor(bias_over, levels = unique(bias_over))
  } else {
    bias_over <- bias_over %||% factor(meta_data[,split_bias_over], levels = unique(meta_data[,split_bias_over]))
  }
  
  # Only keep data that matches bias_over
  your_data <- your_data[,(meta_data[,split_bias_over] %in% bias_over)]
  meta_data <- meta_data[(meta_data[,split_bias_over] %in% bias_over),]
  
  # Order data and metadata
  data_ordered <- your_data[,order(bias_over,meta_data[,split_bias_on])]
  meta_ordered <- meta_data[order(bias_over,meta_data[,split_bias_on]),]
  
  #ensure that filtering results in a subset of samples that is identified by a unique element in plot_over
  if(length(unique(bias_over))*2 != ncol(data_ordered)){
    stop("There should be 2 selections for every unique value of bias_over")
  }

  your_data <- data_ordered[rowSums(data_ordered) > 0,]

  your_data_list <- lapply(2*(1:(ncol(your_data)/2)), function(i){
    temp <- as.data.frame(as.matrix(your_data[,c(i-1,i)])) #subset the data
    temp <- temp[rowSums(temp) > 0,] # only keep barcodes with non-zero values for this pair
    temp <- as.data.frame(prop.table(as.matrix(temp+1), margin = 2)) # Add 1 to all values
    colnames(temp) <- c("F", "S")
    temp$added_prop <- rowSums(temp)
    temp$bias <- temp[,1]/temp[,2] # Calculate bias
    temp$log_bias <- log2(temp$bias)
    temp$TP <- unique(bias_over)[i/2]
    return(temp)
  })
  your_data <- do.call(rbind, your_data_list)

  # Weighted ridge plot
  if (weighted == T){
    as.data.frame(your_data) %>% group_by(TP) %>%
      do(ggplot2:::compute_density(.$log_bias, .$added_prop)) %>%
      dplyr::rename(log_bias = x) -> your_data_densities

    p <- ggplot2::ggplot(your_data_densities, ggplot2::aes(x = log_bias, y = TP, fill = TP, group = TP, height = density)) +
      ggridges::geom_density_ridges(stat="identity", scale = scale)+
      ggplot2::scale_x_continuous(name = paste0("log bias: log2(", bias_1, "/", bias_2, ")")) +
      ggplot2::theme(text = ggplot2::element_text(size = text_size),
                     panel.grid.major.x = ggplot2::element_line(colour = "grey"),
                     panel.grid.major.y = ggplot2::element_blank(),
                     panel.background = ggplot2::element_rect(fill = "white"),
                     axis.text.x = ggplot2::element_text(vjust = 0.5, hjust = 1, angle = 90),
                     legend.position = "none") + ggplot2::scale_y_discrete(name = split_bias_over)
  } else {

  # Normal ridge plot
    
    # Method using ggplot2::compute_diversity - leads to slightly weird scaling
    # as.data.frame(your_data) %>% group_by(TP) %>%
    #   do(ggplot2:::compute_density(.$log_bias, rep(1,times = length(.$log_bias)))) %>%
    #   dplyr::rename(log_bias = x) -> your_data_densities
    
    # p <- ggplot2::ggplot(your_data_densities, ggplot2::aes(x = log_bias, y = TP, fill = TP, group = TP, height = density)) +
    #   ggridges::geom_density_ridges(stat="identity", scale = scale) +
    p <- ggplot2::ggplot(your_data[order(your_data$added_prop),], ggplot2::aes(x = log_bias, y = TP, fill = TP, group = TP)) +
    ggridges::geom_density_ridges(scale = scale) +
    ggplot2::scale_x_continuous(name = paste0("log bias: log2(", bias_1, "/", bias_2, ")")) +
    ggplot2::theme(text = ggplot2::element_text(size = text_size),
                panel.grid.major.x = ggplot2::element_line(colour = "grey"),
                panel.grid.major.y = ggplot2::element_blank(),
                panel.background = ggplot2::element_rect(fill = "white"),
                axis.text.x = ggplot2::element_text(vjust = 0.5, hjust = 1, angle = 90),
                legend.position = "none") + ggplot2::scale_y_discrete(name = split_bias_over) +
                ggplot2::coord_cartesian(clip = "off")
  }

p

}
