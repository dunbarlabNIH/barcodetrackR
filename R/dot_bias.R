#' Dot Bias
#'
#' Given a data frame organized with two alternating cell lineages, gives bias jittered dot plots showing percent total contribution to both lineages over time.
#'@param your_data Your data frame.
#'@param text_size The size of the text in the plot.
#'@param bias_form One of "log2", "fractional", "fold" specifying how to calculate bias.
#'@return Bias plot for two lineages over time.
#'@examples
#'clonal_bias(datum[,c("T 1m", "B 1m", "T 3m" "B 3m")])
#'@export
dot_bias <- function(your_data, text_size = 20, bias_form = "fold")
{
  if (ncol(your_data)%%2 != 0) {
    stop("Data frame must be divisible by 2.")
  }
  your_data <- your_data[rowSums(your_data) > 0,]
  your_data <- as.data.frame(prop.table(as.matrix(your_data), margin = 2))
  your_data_list <- lapply(2 * (1:(ncol(your_data)/2)), function(i) {
    temp <- as.data.frame(prop.table(as.matrix(your_data[,c(i - 1, i)]), margin = 2))
    colnames(temp) <- c("F", "S")
    temp <- temp[rowSums(temp)>0,]
    temp$added_prop <- rowSums(temp)
    if(bias_form == "fold"){
      temp$bias <- ifelse(temp$F >= temp$S, temp$F/temp$S, -temp$S/temp$F)
      temp$bias[is.nan(temp$bias)] <- 1
      temp$bias <- ifelse(temp$bias > 10, 10, temp$bias)
      temp$bias <- ifelse(temp$bias < -10, -10, temp$bias)
      temp$bias <- ifelse(temp$bias < 0, temp$bias + 2, temp$bias)
      temp$TP <- i/2
    }
    return(temp)
  })
  your_data <- do.call(rbind, your_data_list)
  if(bias_form == "fold"){
    format_back = function(x){
      ifelse(x < 1 , (x * -1)+2, x)
    }
    breakies = seq(from = -9, to = 10)
    y_axis_label = "Fold Change"
  }
  ggplot2::ggplot(your_data, ggplot2::aes(x = TP, y = bias))+
    ggplot2::geom_jitter(ggplot2::aes(size = added_prop), color = "blue", fill = "blue", alpha = 0.1) +
    ggplot2::scale_size(range = c(1,20), guide = FALSE)+
    ggplot2::scale_y_continuous(labels = format_back, breaks = breakies)+
    ggplot2::scale_x_continuous(breaks = 1:length(unique(your_data$TP)))+
    ggplot2::ylab(y_axis_label)+
    ggplot2::theme(text = ggplot2::element_text(size = text_size),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(colour = "grey"),
                   panel.background = ggplot2::element_rect(fill = "white", colour = "black"),
                   panel.margin = ggplot2::unit(2,"lines"),
                   axis.text.x = ggplot2::element_text(vjust = 0.5, hjust = 1, angle = 90))

}
