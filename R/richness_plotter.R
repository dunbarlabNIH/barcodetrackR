#' Richness Plotter
#'
#' Plots (using ggplot2) a line plot that tracks the barcode counts of cell types across time points.
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param months Numeric vector. Corresponds to the months of your_data.
#'@param celltypes Character vector. Corresponds to the celltypes of your_data.
#'@param combine Logical. Whether to combine samples from a single month together.
#'@param thresh Numeric. Universal threshold applied to your_data (anything below this threshold will be set to 0)
#'@param point_size Numeric. Size of points.
#'@param line_size Numeric. Size of lines.
#'@param richness_type Character. One of "unique", "cumulative", or "new".
#'@param y_lower Numeric. Lower bound of y axis.
#'@param y_upper Numeric. Upper bound of y axis.
#'@return Outputs plot of diversity tracked over time.
#'@note Data must be organized in columns by month, and then by cell type. For example, when tracking B and T cells
#'over three months, your columns should be in this order: "1m T", "1m B", "2m T", "2m B", "3m T", "3m B". Then
#'celltypes should be c("T", "B"), and months should be c(1,2,3).
#'@examples
#'richness_plot(your_data = zh33_first3m_T_B_NK,
#'                months = c(1,2,3), celltypes = c("T", "B", "NK"),
#'                thresh = 1000, richness_type = "cumulative")
#'
#'@export

richness_plot <- function(your_data, months = c(1:your_data), celltypes, combine = FALSE, thresh = 0, point_size = 5,
                          line_size = 3, richness_type = "unique", y_lower = 0, y_upper = 2000, show_table = FALSE){
  your_data[your_data < thresh] <- 0
  if(ncol(your_data) %% (length(months)*length(celltypes)) != 0){
    stop("Error 2: Enter a correct number of months and correct number of celltypes")
  }
  your_data <- as.data.frame(t(your_data))
  if(combine == TRUE){
    your_data <- split(as.data.frame(your_data), rep(1:length(months), each = length(celltypes)))
    your_data <- lapply(your_data, t)
    your_data <- lapply(your_data, rowSums)
    your_data <- do.call(cbind, your_data)
    your_data <- barcodetrackR::barcodecount(your_data, months = months, count = richness_type, threshold = 0)
    your_data <- as.data.frame(your_data)
    rownames(your_data) <- 1:nrow(your_data)
    colnames(your_data) <- "BARCODES"
    your_data$MONTH <- months
    your_data$CELLTYPE <- "ALL"

  } else {
    your_data <- split(your_data, celltypes)
    your_data <- lapply(your_data, t)
    your_data <- lapply(your_data, barcodetrackR::barcodecount, months = months, count = richness_type, threshold = 0)
    your_data <- lapply(your_data, as.data.frame)
    your_data <- do.call(cbind, your_data)
    colnames(your_data) <- celltypes
    your_data$MONTH <- months
    your_data <- reshape2::melt(your_data, id.vars = "MONTH")
    colnames(your_data)[2:3] <- c("CELLTYPE", "BARCODES")
    your_data <- your_data[,c(3,1,2)]
  }

  if(show_table == TRUE){
    return(your_data)
  } else {
    ggplot2::ggplot(your_data, ggplot2::aes(x = MONTH, y = BARCODES, group = CELLTYPE, colour = CELLTYPE))+
      ggplot2::geom_line(size = line_size)+
      ggplot2::geom_point(size = point_size)+
      ggplot2::ylab(paste0(richness_type, " Barcodes"))+
      ggplot2::xlab("Month")+
      ggplot2::coord_cartesian(ylim = c(y_lower, y_upper))

  }


}
