#' Ridge plot
#'
#' Given a data frame organized with two alternating cell lineages, gives ridge plots showing percent total contribution to both lineages.
# '
#'@param your_data Your data frame.
#'@param text_size The size of the text in the plot.
#'@return Bias plot for two lineages over time.
#'@examples
#'ridge_plot(datum[,c("T 1m", "B 1m", "T 3m" "B 3m")])
#'@export
ridge_plot <- function(your_data, text_size = 20){
  if(ncol(your_data)%% 2 != 0){
    stop("Data frame must be divisible by 2.")
  }
  your_data <- your_data[rowSums(your_data) > 0,]
 # your_data <- as.data.frame(prop.table(as.matrix(your_data), margin = 2)) # Do not need this step
  your_data_list <- lapply(2*(1:(ncol(your_data)/2)), function(i){
    temp <- as.data.frame(as.matrix(your_data[,c(i-1,i)])) #subset the data
    temp <- temp[rowSums(temp) > 0,] # only keep barcodes with non-zero values for this pair
    temp <- as.data.frame(prop.table(as.matrix(temp+1), margin = 2)) # Add 1 to all values
    colnames(temp) <- c("F", "S")
    temp$added_prop <- rowSums(temp)
    temp$bias <- temp[,1]/temp[,2] # Calculate bias
    temp$log_bias <- log2(temp$bias)
    # temp$cuts <- cut(temp$log_bias, breaks = c(-Inf, -2, -1.5, -1, -0.5, 0.5, 1,1.5,  2, Inf), include.lowest = TRUE)
    temp$TP <- i/2
    return(temp)
  })
  your_data <- do.call(rbind, your_data_list)
  ggplot2::ggplot(your_data[order(your_data$added_prop),], ggplot2::aes(x = log_bias, y = TP, fill = TP, group = TP)) +
    ggridges::geom_density_ridges() + 
    ggplot2::theme(text = ggplot2::element_text(size = text_size),
                panel.grid.major.x = ggplot2::element_line(colour = "grey"),
                panel.grid.major.y = ggplot2::element_blank(),
                panel.background = ggplot2::element_rect(fill = "white"),
                axis.text.x = ggplot2::element_text(vjust = 0.5, hjust = 1, angle = 90),
                legend.position = "none")
  
}
