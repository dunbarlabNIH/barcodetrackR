#' Diversity Plotter
#'
#' Plots (using ggplot2) a line plot that tracks the diversity of a lineage of cells across time points.
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param months Numeric vector. Corresponds to the months of your_data.
#'@param celltypes Character vector. Corresponds to the celltypes of your_data.
#'@param thresh Numeric. Universal threshold applied to your_data (anything below this threshold will be set to 0)
#'@param point_size Numeric. Size of points.
#'@param line_size Numeric. Size of lines.
#'@param index_type Character. One of "herfindahl-hirschman", "shannon", "simpson", or "invsimpson".
#'@param y_lower Numeric. Lower bound of y axis.
#'@param y_upper Numeric. Upper bound of y axis.
#'@return Outputs plot of diversity tracked over time.
#'@note Data must be organized in columns by month, and then by cell type. For example, when tracking B and T cells
#'over three months, your columns should be in this order: "1m T", "1m B", "2m T", "2m B", "3m T", "3m B". Then
#'celltypes should be c("T", "B"), and months should be c(1,2,3).
#'@examples
#'diversity_plot(your_data = zh33_first3m_T_B_NK,
#'                months = c(1,2,3), celltypes = c("T", "B", "NK"),
#'                thresh = 1000, index_type = "herfindahl-hirschman")
#'
#'@export
diversity_plot <- function(your_data, months, celltypes, thresh = 0, point_size = 5, line_size = 3, index_type = "shannon", y_lower = 0, y_upper = 9.5) {
  bycelltype <- rep(celltypes, length(months))
  bymonth <- rep(months, each = length(celltypes))
  your_data[your_data < thresh] <- 0
  herfindahl = FALSE
  if(index_type == "herfindahl-hirschman"){
    herfindahl = TRUE
    index_type = "simpson"
  }
  your_data <- vegan::diversity(your_data, index = index_type, MARGIN = 2)
  if(herfindahl){
    your_data <- (1-your_data)
  }

  your_data <- reshape2::melt(as.matrix(your_data))
  your_data$Var2 <- NULL
  your_data$celltype <- bycelltype
  your_data$month <- bymonth

  ggplot2::ggplot(your_data, aes(x = month, y = value, group = celltype, colour = celltype))+
    ggplot2::geom_line(size = line_size)+
    ggplot2::geom_point(size = point_size)+
    ggplot2::ylab(paste0(index_type, " diversity"))+
    ggplot2::xlab("Month")+
    ggplot2::coord_cartesian(ylim = c(y_lower, y_upper))

}
