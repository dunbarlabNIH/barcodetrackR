#' Bias dot plot
#'
#' Given a summarized experiment, gives dot plot of log biases for 2 cell types
# '
#'@param your_SE Your SummarizedExperiment of barcode data and associated metadata
#'@param cell_var The column of metadata corresponding to cell types
#'@param cell_1 The first cell type to be compared
#'@param cell_2 The second cell type to be compared
#'@param filter_by The column of metadata to filter by. e.g. Timepoint
#'@param filter_selection The value of filter_by you want to create histogram for e.g. 2 months
#'@param text_size The size of the text in the plot.
#'@param alpha The transparency of the dots. Lower is more translucent, 1 is completely opaque
#'@return Histogram of log bias for two lineages over time.
#'@examples
#'bias_dot(your_SE = SE, cell_var = "Cell_type", cell_1 = "B", cell_2 = "T", filter_by = "Timepoint", filter_selection = "3m")
#'@export
bias_dot <- function(your_SE, cell_var, cell_1, cell_2, filter_by, filter_selection, text_size = 20, alpha = 0.5){
  # Load data
  your_data <- SummarizedExperiment::assays(your_SE)$counts
  meta_data <- SummarizedExperiment::colData(your_SE)
  # Only keep data that matches filter
  your_data <- your_data[,meta_data[,filter_by] == filter_selection]
  meta_data <- meta_data[meta_data[,filter_by] == filter_selection,]
  # Only keep desired cell types
  your_data <- your_data[,meta_data[,cell_var] == cell_1 | meta_data[,cell_var] == cell_2]
  meta_data <- meta_data[meta_data[,cell_var] == cell_1 | meta_data[,cell_var] == cell_2,]
  # Make the cell_var an ordered factor
  meta_data[,cell_var] <- factor(meta_data[,cell_var], levels = c(cell_1,cell_2))
  # Order data correctly
  your_data <- your_data[,order(meta_data[,cell_var])]
  meta_data <- meta_data[order(meta_data[,cell_var]),]
  
  # Basic error handling
  if(ncol(your_data)%% 2 != 0){
    stop("Data frame must be divisible by 2.")
  }
  your_data <- your_data[rowSums(your_data) > 0,] # Remove rows of data where both cell types have 0 count
  # your_data <- as.data.frame(prop.table(as.matrix(your_data), margin = 2))
  your_data_list <- lapply(2*(1:(ncol(your_data)/2)), function(i){
    temp <- as.data.frame(as.matrix(your_data[,c(i-1,i)])) #subset data
    # temp <- temp[rowSums(temp) > 0,] # Remove b.c.s that are 0 in both pairs
    temp <- as.data.frame(prop.table(as.matrix(temp+1), margin = 2)) # Add 1 to all values
    colnames(temp) <- c("F", "S")
    temp$added_prop <- rowSums(temp)
    temp$bias <- temp[,1]/temp[,2]
    temp$log_bias <- log2(temp$bias)
    temp$cuts <- cut(temp$log_bias, breaks = seq(from = floor(min(as.numeric(temp$log_bias))),to = ceiling(max(as.numeric(temp$log_bias))),by = (ceiling(max(as.numeric(temp$log_bias)))-floor(min(as.numeric((temp$log_bias)))))*1/10), include.lowest = TRUE)
    temp$TP <- i/2
    return(temp)
  })
  your_data <- do.call(rbind, your_data_list)

  # return(floor(min(your_data$log_bias)))
  ggplot2::ggplot(your_data, ggplot2::aes(x = TP, y = log_bias))+
    ggplot2::geom_jitter(ggplot2::aes(size = added_prop), color = "blue", fill = "blue", alpha = alpha) +
  #  ggplot2::scale_size(range = c(floor(min(your_data$log_bias)),ceiling(max(your_data$log_bias))), guide = FALSE)+
    # ggplot2::scale_y_continuous(breaks = your_data$cuts)+
    # ggplot2::scale_x_continuous(breaks = 1:length(unique(your_data$TP)))+
    ggplot2::ylab("log_bias")+
    ggplot2::theme(text = ggplot2::element_text(size = text_size),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(colour = "grey"),
                   panel.background = ggplot2::element_rect(fill = "white", colour = "black"),
                   panel.spacing = ggplot2::unit(2,"lines"),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank()) +
    ggplot2::labs(size = "Proportion", x = "", title = paste(filter_by,"=",filter_selection))+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))+
    ggplot2::coord_cartesian(clip = "off")+
    ggplot2::annotate("text",x = Inf, y = Inf, label = cell_1, size = 8,hjust = -1, vjust = 1)+
    ggplot2::annotate("text",x = Inf, y = -Inf, label = cell_2, size = 8, hjust = -1, vjust = 0)


}


