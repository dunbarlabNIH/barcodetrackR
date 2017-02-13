#'@title Ternary Plot
#'
#'@description Creates a ternary plot showing the bias of clones towards one of three axes.
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param show_arrows Logical. Display arrows or not.
#'@param show_ticks Logical. Display tick marks or not.
#'@param density_mode Logical. Uses a kernel density estimation to view concentrations of clones.
#'@return Displays a ternary plot in the current plot window.
#'@examples
#'ternary_plot(your_data = zh33, dot_size = 1000, density_mode = TRUE)
#'@export
ternary_plot <- function(your_data, show_arrows = TRUE, show_ticks = TRUE, density_mode = FALSE){
  if (ncol(your_data) != 3){
    stop("You must include 3 columns of data.")
  }

  label_names <- paste0(colnames(your_data))
  colnames(your_data) <- c("X1", "X2", "X3")
  your_data <- your_data[rowSums(your_data) != 0,]
  your_data <- as.data.frame(prop.table(as.matrix(your_data), margin = 2))
  ABUNDANCE <- rowSums(your_data)
  your_data <- as.data.frame(prop.table(as.matrix(your_data), margin = 1))

  if(density_mode){
    your_data$ABUNDANCE <- 1
    your_data$MYCOLORS <- 'black'
  } else {
    your_data$ABUNDANCE <- ABUNDANCE
    your_data$MYCOLORS <- rainbow(nrow(your_data),s = 1, v = 1, start = 0, end = 0.75, alpha = 1)
  }

  your_data <- your_data[order(-your_data$ABUNDANCE),]


  g <- ggtern::ggtern(your_data, ggtern::aes(X1, X2, X3))+
    ggtern::coord_tern(expand = TRUE)+
    ggtern::limit_tern(breaks = seq(0,1,by=.20), labels = paste0(c(0,20,40,60,80,100), "%"))+
    ggplot2::labs(x = label_names[1], y = label_names[2], z = label_names[3])+
    ggtern::theme(
      tern.axis.ticks.length.major = ggplot2::unit(30, units = "points"),
      tern.axis.text =ggplot2::element_text(vjust=0.3, colour="black", size = 15)

    )
    if (show_arrows == TRUE)
      g <- g + ggtern::theme_arrowlong()
  if (show_ticks == FALSE)
    g <- g + ggtern::theme_hidelabels()+ggtern::theme_hideticks()

  if (density_mode == TRUE) {
    g <- g + ggtern::stat_density_tern(geom='polygon', ggtern::aes(fill=..level..), base="identity", colour = 'black') +
      ggplot2::scale_fill_gradient(low='green',high='red', guide = FALSE)+
      ggtern::geom_mask()+
      ggplot2::geom_point(ggplot2::aes(size = ABUNDANCE), color = "black", show.legend = FALSE)
  } else {
    g <- g + ggtern::geom_mask()+
      ggplot2::geom_point(ggplot2::aes(size = ABUNDANCE, fill = MYCOLORS), shape = 21)+
      ggplot2::scale_fill_discrete(guide = FALSE)
  }

  par(mar=c(5,6,40,2))
  g


}
