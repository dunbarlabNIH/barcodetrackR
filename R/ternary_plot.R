ternary_plot <- function(your_data, colors = "random", dot_size = 2000, scale_abundance = TRUE,
                         show_numeric_scale = FALSE, corner_labels = 20, show_arrows = TRUE,
                         show_breaks = TRUE, density_mode = FALSE){
  if (ncol(your_data) != 3){
    stop("You must include 3 columns of data.")
  }

  your_data <- your_data[rowSums(your_data) != 0,]
  ABUNDANCE <- rowSums(your_data)[rowSums(your_data) != 0]
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
    point_size = ABUNDANCE * dot_size
    point_color = COLORS
  }


  g <- ggtern::ggtern(your_data, ggtern::aes(X1, X2, X3))+
    ggplot2::geom_point(size = point_size, fill = point_color, shape = 21)+
    ggplot2::labs(x = label_names[1], y = label_names[2], z = label_names[3])
  if (density_mode == TRUE) {
    g <- g + ggtern::stat_density_tern(geom='polygon', ggtern::aes(fill=..level..), base="identity", colour = 'black') +
      ggplot2::scale_fill_gradient(low='green',high='red', guide = FALSE) + ggtern::theme_bw(base_size = corner_labels) + ggtern::theme_nogrid()
  } else {
    g <- g + ggplot2::theme_minimal(base_size = corner_labels)
  }
  if (show_arrows == TRUE)
    g <- g + ggtern::theme_arrowlong()

  g+ggtern::theme_hidelabels()+ggtern::theme_hideticks()


}
