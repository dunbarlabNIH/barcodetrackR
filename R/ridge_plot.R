#' Ridge plot
#'
#' Given a summarized experiment, gives ridge plots showing percent total contribution to both lineages.
# '
#'@param your_SE Your SummarizedExperiment of barcode data and associated metadata
#'@param cell_var The column of metadata corresponding to cell types
#'@param cell_1 The first cell type to be compared. Will be on the RIGHT side of the ridge plot
#'@param cell_2 The second cell type to be compared. Will be on the LEFT side of the ridge plot
#'@param plot_by The column of metadata to plot by. If numeric, y axis will be in increasing order. If categorical, it will follow order of metadata
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
#'ridge_plot(your_se = SE, cell_var = "selection_type", selection_1 = "B", selection_2 = "T", plot_by = "Timepoint")
#'@export
ridge_plot <- function(your_SE,
                       selection_var,
                       selection_1,
                       selection_2,
                       plot_by,
                       weighted = F,
                       text_size = 16,
                       scale = 1){
  
  # Load data
  your_data <- SummarizedExperiment::assays(your_SE)$counts
  meta_data <- SummarizedExperiment::colData(your_SE)

  # Basic error handling
  coldata_names <- colnames(SummarizedExperiment::colData(your_SE))
  if(any(! c(selection_var) %in% coldata_names)){
    stop("selection_var must match a column name in colData(your_SE)")
  }
  if(! selection_1 %in% unique(SummarizedExperiment::colData(your_SE)[[selection_var]])){
    stop("selection_1 is not in selection_var")
  }
  if(! selection_2 %in% unique(SummarizedExperiment::colData(your_SE)[[selection_var]])){
    stop("selection_2 is not in selection_var")
  }
  
  # Only keep data that matches filter
  your_data <- your_data[,meta_data[,selection_var] == selection_1 | meta_data[,selection_var] == selection_2]
  meta_data <- meta_data[meta_data[,selection_var] == selection_1 | meta_data[,selection_var] == selection_2,]
  
  # Make the selection_var an ordered factor
  meta_data[,selection_var] <- factor(meta_data[,selection_var], levels = c(selection_1,selection_2))
  
  # error handling
  if(any(! c(plot_by) %in% coldata_names)){
    stop("plot_by must match a column name in colData(your_SE)")
  }
  
  # Order data columns based on meta data to plot by
  if (class(meta_data[,plot_by]) == "numeric"){
    meta_ordered <- meta_data[order(meta_data[,plot_by],meta_data[,selection_var]),]
    # Change to factor
    meta_ordered[,plot_by] <- factor(meta_ordered[,plot_by])
    data_ordered <- your_data[,order(meta_data[,plot_by],meta_data[,selection_var])]
  } else {
    # Order meta_data based on the order of the data not alphabetically
    meta_data[,plot_by] <- factor(meta_data[,plot_by], levels = unique(meta_data[,plot_by]))
    meta_ordered <- meta_data[order(meta_data[,plot_by],meta_data[,selection_var]),]
    data_ordered <- your_data[,order(meta_data[,plot_by],meta_data[,selection_var])]
  }
  
  #ensure that filtering results in a subset of samples that is identified by a unique element in plot_over
  if(length(unique(meta_ordered[,plot_by]))*2 != ncol(data_ordered)){
    stop("There should be 2 selections for every unique value of plot_by")
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
    temp$TP <- unique(meta_ordered[,plot_by])[i/2]
    return(temp)
  })
  your_data <- do.call(rbind, your_data_list)

  # Weighted ridge plot
  if (weighted == T){
    as.data.frame(your_data) %>% group_by(TP) %>%
      do(ggplot2:::compute_density(.$log_bias, .$added_prop)) %>%
      dplyr::rename(log_bias = x) -> your_data_densities

    p <- ggplot2::ggplot(your_data_densities, ggplot2::aes(x = log_bias, y = TP, fill = TP, group = TP, height = density)) +
      ggridges::geom_density_ridges(stat="identity", scale = scale) +
      ggplot2::theme(text = ggplot2::element_text(size = text_size),
                     panel.grid.major.x = ggplot2::element_line(colour = "grey"),
                     panel.grid.major.y = ggplot2::element_blank(),
                     panel.background = ggplot2::element_rect(fill = "white"),
                     axis.text.x = ggplot2::element_text(vjust = 0.5, hjust = 1, angle = 90),
                     legend.position = "none") + ggplot2::scale_y_discrete(name = plot_by)
                     # ggplot2::annotate("text",x = Inf, y = 0.6, label = selection_1, size = 8,hjust = 5)+
                     # ggplot2::annotate("text",x = -Inf, y = 0.6, label = selection_2, size = 8, hjust = -5)
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
    ggplot2::theme(text = ggplot2::element_text(size = text_size),
                panel.grid.major.x = ggplot2::element_line(colour = "grey"),
                panel.grid.major.y = ggplot2::element_blank(),
                panel.background = ggplot2::element_rect(fill = "white"),
                axis.text.x = ggplot2::element_text(vjust = 0.5, hjust = 1, angle = 90),
                legend.position = "none") + ggplot2::scale_y_discrete(name = plot_by) +
                ggplot2::coord_cartesian(clip = "off")
                # ggplot2::annotate("text",x = Inf, y = -Inf, label = selection_1, size = 8,hjust = 5, vjust = -1)+
                # ggplot2::annotate("text",x = -Inf, y = -Inf, label = selection_2, size = 8, hjust = -7, vjust = -1)
                # # The hjust above needs to be able to scale to different x limits
  }

p

}
