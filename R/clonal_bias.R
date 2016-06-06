#' Clonal Bias
#'
#' Given a data frame organized with two alternating cell lineages, gives bias plots showing percent total contribution to both lineages.
#'@param your_data Your data frame.
#'@param text_size The size of the text in the plot.
#'@return Bias plot for two lineages over time.
#'@examples
#'clonal_bias(datum[,c("T 1m", "B 1m", "T 3m" "B 3m")])
#'@export
clonal_bias <- function(your_data, text_size = 20){
  if(ncol(your_data)%% 2 != 0){
    stop("Data frame must be divisible by 2.")
  }
  your_data <- your_data[rowSums(your_data) > 0,]
  your_data <- as.data.frame(prop.table(as.matrix(your_data), margin = 2))
  your_data_list <- lapply(2*(1:(ncol(your_data)/2)), function(i){
    temp <- as.data.frame(prop.table(as.matrix(your_data[,c(i-1,i)]), margin = 2))
    temp <- temp[rowSums(temp) > 0,]
    colnames(temp) <- c("F", "S")
    temp$added_prop <- rowSums(temp)
    temp$bias <- temp[,1]/temp[,2]
    temp$log_bias <- log2(temp$bias)
    temp$cuts <- cut(temp$log_bias, breaks = c(-Inf, -2, -1.5, -1, -0.5, 0.5, 1,1.5,  2, Inf), include.lowest = TRUE)
    temp$TP <- i/2
    return(temp)
  })
  your_data <- do.call(rbind, your_data_list)
  ggplot2::ggplot(your_data[order(your_data$added_prop),], ggplot2::aes(x = cuts, y = added_prop, color = cuts))+
    ggplot2::facet_wrap(~ TP)+
    ggplot2::scale_color_manual(values = colorRampPalette(c("navyblue", "springgreen", "navyblue"))(length(unique(your_data$cuts))))+
    ggplot2::geom_bar(stat = "identity", fill = "white")+
    ggplot2::scale_x_discrete()+
    ggplot2::theme(text = ggplot2::element_text(size = text_size),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(colour = "grey"),
                   panel.background = ggplot2::element_rect(fill = "white", colour = "black"),
                   panel.margin = ggplot2::unit(2, "lines"),
                   axis.text.x = ggplot2::element_text(vjust = 0.5, hjust = 1, angle = 90))




}


