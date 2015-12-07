#' Richness Plotter
#'
#' Plots (using ggplot2) a line plot that tracks the barcode counts of cell types across time points.
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param months Numeric vector. Corresponds to the months of your_data.
#'@param celltypes Character vector. Corresponds to the celltypes of your_data.
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
#'

richness_plot <- function(your_data, months, celltypes, thresh = 0, point_size = 5, line_size = 3, richness_type = "unique", y_lower = 0, y_upper = 2000){
  your_data[your_data < thresh] <- 0
  dfr <- as.data.frame(t(your_data))
  dfr <- split(dfr, celltypes)
  dfr <- lapply(dfr, t)
  dfr <- lapply(dfr, barcodetrackR::barcodecount, months = months, count = richness_type, threshold = 0)
  newnames <- rep(months, length(celltypes))
  newcellorder <- names(dfr)
  newnames <- paste0(months, "m ", rep(names(dfr), each = length(months)))
  dfr <- unlist(dfr)
  names(dfr) <- newnames
  dfr <- reshape2::melt(dfr)
  dfr$celltype <- rep(newcellorder, each = length(months))
  dfr$month <- rep(months, length(celltypes))
  ggplot2::ggplot(dfr, ggplot2::aes(x = month, y = value, group = celltype, colour = celltype))+
    ggplot2::geom_line(size = line_size)+
    ggplot2::geom_point(size = point_size)+
    ggplot2::ylab(paste0(richness_type, " Barcodes"))+
    ggplot2::xlab("Month")+
    ggplot2::coord_cartesian(ylim = c(y_lower, y_upper))
  }
