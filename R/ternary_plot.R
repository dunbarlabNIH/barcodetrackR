#' #'@title Ternary Plot
#' #'
#' #'@description Creates a ternary plot showing the bias of clones towards one of three axes.
#' #'
#' #'@param your_SE A SummarizedExperiment.
#' #'@param show_arrows Logical. Display arrows or not.
#' #'@param show_ticks Logical. Display tick marks or not.
#' #'@param density_mode Logical. Uses a kernel density estimation to view concentrations of clones.
#' #'@return Displays a ternary plot in the current plot window.
#' #'@examples
#' #'ternary_plot(your_SE = zh33, density_mode = TRUE)
#' ternary_plot <- function(your_SE, show_arrows = FALSE, show_ticks = FALSE, density_mode = TRUE){
#' 
#'   # Load data
#'   your_data <- SummarizedExperiment::assays(your_SE)[["percentages"]]
#'   meta_data <- SummarizedExperiment::colData(your_SE)
#'   your_data <- your_data[rowSums(your_data) > 0,]
#'   label_names <- colnames(your_data)
#'   if (ncol(your_data) != 3){
#'     stop("You must include exactly 3 columns of data.")
#'   }
#' 
#'   your_data %>%
#'     magrittr::set_colnames(c("X1", "X2", "X3")) %>%
#'     dplyr::mutate(cumulative = sum(X1, X2, X3)) %>%
#'     dplyr::mutate(average_cumulative = cumulative/3) %>%
#'     dplyr::mutate(avg_1 = X1/cumulative, avg_2 = X2/cumulative, avg_3 = X3/cumulative) -> plotting_data
#' 
#'   # if(density_mode){
#'   #   your_data$ABUNDANCE <- 1
#'   #   your_data$MYCOLORS <- 'black'
#'   # } else {
#'   #   your_data$ABUNDANCE <- ABUNDANCE
#'   #   your_data$MYCOLORS <- rainbow(nrow(your_data),s = 1, v = 1, start = 0, end = 0.75, alpha = 1)
#'   # }
#'   # your_data <- your_data[order(-your_data$ABUNDANCE),]
#' 
#' 
#'   g <- ggtern::ggtern(data = plotting_data, ggtern::aes(x = X1, y = X2, z = X3))+
#'     ggtern::coord_tern(expand = TRUE)+
#'     ggtern::limit_tern(breaks = seq(0,1,by=.20), labels = paste0(c(0,20,40,60,80,100), "%"))+
#'     ggplot2::labs(x = label_names[1], y = label_names[2], z = label_names[3])+
#'     ggtern::theme(
#'       tern.axis.ticks.length.major = ggplot2::unit(30, units = "points"),
#'       tern.axis.text =ggplot2::element_text(vjust=0.3, colour="black", size = 15),
#'       tern.axis.arrow = ggplot2::element_line(color = "black", size = 5,),
#'       tern.axis.title.L = ggplot2::element_text(size = 8, angle = -45),
#'       tern.axis.title.R = ggplot2::element_text(size = 8, angle = 45, hjust = 0.5),
#'       tern.axis.title.T = ggplot2::element_text(size = 8),
#'     )
#' 
#' 
#'   # if (show_arrows == TRUE)
#'   #   g <- g + ggtern::theme_arrowlong()
#'   # if (show_ticks == FALSE)
#'   #   g <- g + ggtern::theme_hidelabels()+ggtern::theme_hideticks()
#'   # if (density_mode == TRUE) {
#'   #   g <- g + ggtern::stat_density_tern(geom='polygon', ggtern::aes(fill=..level..), base="identity", colour = 'black') +
#'   #     ggplot2::scale_fill_gradient(low='green', high='red', guide = FALSE)+
#'   #     ggtern::geom_mask() #+
#'   #   #ggplot2::geom_point(ggplot2::aes(size = ABUNDANCE), color = "black", show.legend = FALSE)
#'   # } else {
#'   #   g <- g + ggtern::geom_mask()+
#'   #     ggplot2::geom_point(ggplot2::aes(size = ABUNDANCE, fill = MYCOLORS), shape = 21)+
#'   #     ggplot2::scale_fill_discrete(guide = FALSE)+
#'   #     ggplot2::scale_size_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1),
#'   #                                    limits = c(0,1),
#'   #                                    range = c(0,10),
#'   #                                    labels = c("0%", "20%", "40%", "60%", "80%", "100%"))
#'   # }
#' 
#'   g
#' 
#' 
#' }
