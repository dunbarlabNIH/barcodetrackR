ternary_plot <- function(your_data, colors = "random", dot_size = 2000, scale_abundance = TRUE, show_numeric_scale = FALSE, corner_labels = 1){
  if (ncol(your_data) != 3){
    stop("You must include 3 columns of data.")
  }
  ABUNDANCE <- rowSums(your_data)

  your_data <- your_data[ABUNDANCE != 0,]
  ABUNDANCE <- ABUNDANCE[ABUNDANCE != 0]

  ABUNDANCE <- ABUNDANCE/sum(ABUNDANCE)

  your_data <- prop.table(as.matrix(your_data), margin = 2)
  your_data <- as.data.frame(prop.table(as.matrix(your_data), margin = 1))

  COLORS <- rainbow(nrow(your_data),s = 1, v = 1, start = 0, end = 0.75, alpha = 1)

  label_names <- colnames(your_data)
  colnames(your_data) <- c("X1", "X2", "X3")


  ggtern::ggtern(your_data, ggtern::aes(X1, X2, X3))+
    ggplot2::geom_point(fill=COLORS,shape=21,size=dot_size*ABUNDANCE)+
    ggplot2::theme_minimal()+
    ggplot2::labs(x = label_names[1], y = label_names[2], z = label_names[3], size = 10)
}
