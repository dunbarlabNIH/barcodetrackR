#'@title Ternary Plot
#'
#'@description Creates a ternary plot showing the bias of clones towards one of three axes.
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param dot_size The size of the dots (does not apply to density mode).
#'@param show_arrows Logical. Display arrows or not.
#'@param show_ticks Logical. Display tick marks or not.
#'@param density_mode Logical. Uses a kernel density estimation to view concentrations of clones.
#'@return Displays a ternary plot in the current plot window.
#'@examples
#'ternary_plot(your_data = zh33, dot_size = 1000, density_mode = TRUE)
#'@export
ternary_plot <- function(your_data, dot_size = 2000, show_arrows = TRUE, show_ticks = TRUE, density_mode = FALSE){
  if (ncol(your_data) != 3){
    stop("You must include 3 columns of data.")
  }

  your_data <- your_data[rowSums(your_data) != 0,]
  ABUNDANCE <- rowSums(your_data)
  ABUNDANCE <- ABUNDANCE/sum(ABUNDANCE)
  your_data <- prop.table(as.matrix(your_data), margin = 2)
  your_data <- as.data.frame(prop.table(as.matrix(your_data), margin = 1))

  COLORS <- rainbow(nrow(your_data),s = 1, v = 1, start = 0, end = 0.75, alpha = 1)

  label_names <- colnames(your_data)
  colnames(your_data) <- c("X1", "X2", "X3")



  if(density_mode == TRUE){
    point_size = .5
    point_color = "black"
  } else {
    your_data <- your_data[order(-ABUNDANCE),]
    point_color <- COLORS[order(-ABUNDANCE)]
    point_size <- (ABUNDANCE*dot_size)[order(-ABUNDANCE)]
  }


  g <- ggtern::ggtern(your_data, ggtern::aes(X1, X2, X3))
  if (density_mode == TRUE) {
    g <- g + ggtern::stat_density_tern(geom='polygon', ggtern::aes(fill=..level..), base="identity", colour = 'black') +
      ggplot2::scale_fill_gradient(low='green',high='red', guide = FALSE) + ggtern::theme_bw() + ggtern::theme_nogrid()
  } else {
    g <- g + ggplot2::theme_minimal()
  }
  g <- g + ggplot2::geom_point(size = point_size, fill = point_color, shape = 21)+
    ggplot2::labs(x = label_names[1], y = label_names[2], z = label_names[3])
  if (show_arrows == TRUE)
    g <- g + ggtern::theme_arrowlong()
  if (show_ticks == FALSE)
    g <- g + ggtern::theme_hidelabels()+ggtern::theme_hideticks()
  g


}
