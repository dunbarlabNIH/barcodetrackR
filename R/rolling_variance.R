rolling_variance <- function(your_data, text_size = 20) 
{
  if (ncol(your_data)%%2 != 0) {
    stop("Data frame must be divisible by 2.")
  }
  your_data <- your_data[rowSums(your_data) > 0, ]
  your_data <- as.data.frame(prop.table(as.matrix(your_data), margin = 2))
  your_data_added <- as.data.frame(rowSums(your_data))
  your_data_added$barcode <- rownames(your_data_added)
  colnames(your_data_added)[1] <- "sum"
  your_data_added <- your_data_added[order(your_data_added$sum),]
  your_data_added$ordered <- as.factor(1:nrow(your_data_added))
  print(head(your_data_added))
  your_data_list <- lapply(2 * (1:(ncol(your_data)/2)), function(i) {
    temp <- as.data.frame(prop.table(as.matrix(your_data[,c(i - 1, i)]), margin = 2))
    #     temp <- temp[rowSums(temp) > 0, ]
    colnames(temp) <- c("F", "S")
    temp$added_prop <- rowSums(temp)
    temp$bias <- temp[, 2]/(temp[, 1] + temp[, 2])
    temp$bias[is.nan(temp$bias)] <- 0.5
    # temp$log_bias <- log2(temp$bias)
    temp$TP <- i/2
    temp$barcode <- rownames(temp)
    return(temp)
  })
  your_data <- do.call(rbind, your_data_list)
  print(dim(your_data))
  your_data_bias <- reshape2::dcast(your_data, barcode ~ TP, value.var = "bias")
  z <- t(your_data_bias[2:ncol(your_data_bias)])
  print("h")
  something <- zoo::rollapplyr(z, width = 1:ncol(your_data_bias), FUN = var)
  print("j")
  something <- as.data.frame(t(something))
  something[is.na(something)] <- 0
  colnames(something) <- 1:(ncol(your_data_bias)-1)
  something$barcode <- your_data_bias$barcode
  something <- reshape2::melt(something, id.vars = "barcode")
  colnames(something) <- c("barcode", "TP", "rolling")
  something <- merge(something, your_data, by = c("TP", "barcode"))
  something <- merge(something, your_data_added, by = "barcode")
  colnames(something)[8] <- "allsum"
  # return(something)
  print(head(something))
  ggplot2::ggplot(something, ggplot2::aes(x = TP, y = rolling))+ 
    ggplot2::geom_line(ggplot2::aes(color = allsum, group = ordered))+
    ggplot2::scale_color_gradient(low = "grey", high = "black")+
    ggplot2::theme(text = ggplot2::element_text(size = text_size), 
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(colour = "grey"), 
                   panel.background = ggplot2::element_rect(fill = "white", colour = "black"),
                   panel.margin = ggplot2::unit(2,"lines"),
                   axis.text.x = ggplot2::element_text(vjust = 0.5, hjust = 1, angle = 90))
}