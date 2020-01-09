#' Bias histogram
#'
#' Given a summarized experiment, gives histogram of log biases for 2 cell types
# '
#'@param your_SE Your SummarizedExperiment of barcode data and associated metadata
#'@param split_bias_on The column in `colData(your_SE)` from which `bias_1` and `bias_2` will be chosen
#'@param bias_1 The factor you wish to plot on the left side of the plots.
#'@param bias_2 The factor you wish to plot on the right side of the plots.
#'@param split_bias_over The column in `colData(your_SE)` that you wish to observe the bias split on. Defaults to all.
#'@param bias_over Choice(s) from the column designated in `split_bias_over` that will be used for plotting.
#'@param breaks Numeric. The breaks specified for bins on the x-axis (how biased the clones are towards one factor or the other).
#'@param text_size The size of the text in the plot.
#'@param linesize The linewidth of the stacked bars which represent individual barcodes
#'@return Histogram of log bias for two lineages over time.
#'@examples
#'bias_histogram(your_SE = SE, split_bias_on = "Cell_type", bias_1 = "B", bias_2 = "T", split_bias_over = "Timepoint", bias_over = "3m")
#'@export
bias_histogram <- function(your_SE, split_bias_on, bias_1, bias_2, split_bias_over, bias_over = NULL, breaks = c(10,5,2,1), text_size = 20, linesize = .4){

  # Some basic error checking before running the function
  if(length(breaks) != length(unique(breaks))){
    stop("breaks must be unique")
  }
  coldata_names <- colnames(SummarizedExperiment::colData(your_SE))
  if(any(! c(split_bias_on, split_bias_over) %in% coldata_names)){
    stop("split_bias_on and split_bias_over must both match a column name in colData(your_SE)")
  }
  if(any(! c(bias_1, bias_2) %in% levels(SummarizedExperiment::colData(your_SE)[[split_bias_on]]))){
    stop("bias_1 and bias_2 must both be levels in the colData column specified with split_bias_on")
  }
  bias_over <- bias_over %||% levels(SummarizedExperiment::colData(your_SE)[[split_bias_over]])


  plot_list <- lapply(1:length(bias_over), function(i){
    temp_bias_over <- bias_over[i]
    temp_subset <- your_SE[,(your_SE$split_bias_over %in% temp_bias_over) & (your_SE$split_bias_on %in% c(bias_1, bias_2))]
    print(temp_subset)
  })



  # Load data]
  your_data <- SummarizedExperiment::assays(your_SE)[["percentages"]]
  meta_data <- SummarizedExperiment::colData(your_SE)
  # Only keep data that matches filter
  your_data <- your_data[,meta_data[,split_bias_over] == bias_over]
  meta_data <- meta_data[meta_data[,split_bias_over] == bias_over,]
  # Only keep desired cell types
  your_data <- your_data[,meta_data[,split_bias_on] == bias_1 | meta_data[,split_bias_on] == bias_2]
  meta_data <- meta_data[meta_data[,split_bias_on] == bias_1 | meta_data[,split_bias_on] == bias_2,]
  # Make the split_bias_on an ordered factor
  meta_data[,split_bias_on] <- factor(meta_data[,split_bias_on], levels = c(bias_1,bias_2))
  # Order data correctly
  your_data <- your_data[,order(meta_data[,split_bias_on])]
  meta_data <- meta_data[order(meta_data[,split_bias_on]),]

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
    # xmin <- min(temp$log_bias)
    # xmax <- max(temp$log(bias))
    # my_span <- c(xmin,xmin + (xmax-xmin)*1/12,xmin + (xmax-xmin)*2/12,xmin + (xmax-xmin)*3/12,xmin + (xmax-xmin)*4/12,xmin + (xmax-xmin)*5/12,xmin + (xmax-xmin)*6/12,xmin + (xmax-xmin)*7/12,xmin + (xmax-xmin)*8/12,xmin + (xmax-xmin)*9/12,xmin + (xmax-xmin)*10/12,xmin + (xmax-xmin)*11/12,xmax)
    temp$cuts <- cut(temp$log_bias, breaks = seq(from = floor(min(temp$log_bias)),to = ceiling(max(temp$log_bias)),by = (ceiling(max(temp$log_bias))-floor(min(temp$log_bias)))*1/10), include.lowest = TRUE)
    temp$TP <- i/2
    return(temp)
  })
  your_data <- do.call(rbind, your_data_list)
  ggplot2::ggplot(your_data[order(your_data$added_prop),], ggplot2::aes(x = cuts, y = added_prop, color = cuts))+
    # ggplot2::facet_wrap(~ TP)+
    ggplot2::scale_color_manual(values = colorRampPalette(c("navyblue", "springgreen", "navyblue"))(length(unique(your_data$cuts))))+
    ggplot2::geom_bar(stat = "identity", fill = "white", size = linesize)+
    ggplot2::scale_x_discrete(name = "log_bias")+
    ggplot2::theme(text = ggplot2::element_text(size = text_size),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(colour = "grey"),
                   panel.background = ggplot2::element_rect(fill = "white", colour = "black"),
                   panel.spacing = ggplot2::unit(2, "lines"),
                   axis.text.x = ggplot2::element_text(vjust = 0.5, hjust = 1, angle = 90))+
                   ggplot2::labs(color = "", y = "Proportion",title = paste(split_bias_over,"=",bias_over))+
                   ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))+
                   ggplot2::coord_cartesian(clip = "off")+
                   ggplot2::annotate("text",x = Inf, y = -Inf, label = bias_1, size = 8,hjust = -1, vjust = 2)+
                   ggplot2::annotate("text",x = -Inf, y = -Inf, label = bias_2, size = 8, hjust = 1, vjust = 2)




}


