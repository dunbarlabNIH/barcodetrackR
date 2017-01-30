#'@title Barcode Binary Heat Map
#'
#'@description Creates a binary heatmap showing the presence of new clones from L to R in the dataset.
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param names Character vetor
#'@param threshold A chosen (percentage) threshold that decides if a barcode exists or not.
#'@param label_size The size of the column labels.
#'@param your_title The title.
#'@return Displays a binary heat map in the current plot window.
#'@examples
#'barcode_binary_heatmap(your_data = zh33, percent_threshold = 0.01, label_size = 3)
#'BBHM(your_data = zh33)
#'@export


barcode_binary_heatmap <- function(your_data,
                                   names = colnames(your_data),
                                   your_title = "",
                                   label_size = 1,
                                   percent_threshold = 0) {
  your_data <- as.data.frame(prop.table(as.matrix(your_data),2))
  your_data[is.na(your_data)] <- 0
  if (any(colSums(your_data) == 0)){
    stop("One of your columns contained no data")
  }
  your_data[your_data < percent_threshold] <- 0
  your_data[your_data > 0] <- 1
  your_data <- your_data[rowSums(your_data) > 0,]
  barcode_order <- do.call(order, as.data.frame(your_data))
  plotting_data <- reshape2::melt(as.matrix(your_data))
  colnames(plotting_data) <- c("BARCODE", "SAMPLE", "PRESENCE")
  plotting_data$BARCODE <- factor(plotting_data$BARCODE, levels = rev(rownames(your_data[barcode_order,])))
  plotting_data$PRESENCE <- factor(plotting_data$PRESENCE)
  ggplot2::ggplot(plotting_data, ggplot2::aes(x = SAMPLE, y = BARCODE))+
    ggplot2::geom_tile(ggplot2::aes(fill = PRESENCE))+
    ggplot2::scale_fill_manual(values = c("1" = "grey", "0" = "white"), expand = c(0,0), labels = c("Absent", "Present"))+
    ggplot2::scale_y_discrete(labels = NULL, breaks = NULL, expand = c(0,0))+
    ggplot2::scale_x_discrete(expand = c(0,0), labels = names)+
    ggplot2::ylab(NULL)+
    ggplot2::xlab(NULL)+
    ggplot2::ggtitle(paste0("\n", your_title, "\n"))+
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 20),
      axis.text.x = ggplot2::element_text(angle=90, hjust = 1, vjust = 0.5, size = label_size),
      legend.text = ggplot2::element_text(size =  15, face = 'bold'),
      legend.title = ggplot2::element_text(size =  15),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1)
    )







}
