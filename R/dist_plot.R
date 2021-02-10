#' Pairwise Distance Plot
#'
#' Plots the pairwise distances of the specified assay between each sample-sample pair in the provided SummarizedExperiment.
#'
#'@param your_SE A Summarized Experiment object.
#'@param assay The choice of assay to use for the correlation calculation. Set to "proportions" by default.
#'@param plot_labels Vector of x axis labels. Defaults to colnames(your_SE).
#'@param dist_method Character. Distance OR similarity measure from the `proxy` package. Full list of distance and similarity measures can be found using `summary(proxy::pr_DB)`. Default is "euclidean".
#'Distances will be calculated for distance measures, while similarities will be calculated for similarity measures. Distance OR similarity measure will be calculated using the `assay`
#'specified.
#'@param cluster_tree Logical. Whether to cluster samples and plot a hierarchical tree calculated from the distance or similarity measure used. Default is FALSE.
#'@param your_title Character. The title for the plot.
#'@param grid Logical. Include a grid or not in the resulting plot.
#'@param label_size Numeric. The size of the column labels.
#'@param plot_type Character. One of "color", "circle", or "number".
#'@param no_negatives Logical. Whether to make negative correlations = 0.
#'@param return_table Logical. Whether or not to return table of p-values, confidence intervals, and R values instead of displaying a plot.
#'@param color_pal Character. One of 'Reds', 'Purples', 'Oranges', 'Greys', 'Greens', or 'Blues' that designates the brewer.pal color scale to use.
#'@param number_size Numeric. size of the text label when plot_type is "number".
#'@param point_scale Numeric. The size of the largest point if the plot_type is "circle".
#'@param minkowski_p Numeric. If 'Minkowski' is chosen, the 'p' used to calculate the Minkowski distance.
#'
#'@return Plots pairwise correlation plot for the samples in your_SE.
#'
#'@importFrom rlang %||%
#'@importFrom magrittr %>%
#'
#'@export
#'
#'@examples
#'dist_plot(your_SE = wu_subset, your_title = "Pairwise Euclidean distances between all samples", plot_type = "color")
#"

dist_plot = function(your_SE,
                    assay = "proportions",
                    plot_labels = colnames(your_SE),
                    dist_method = "euclidean",
                    cluster_tree = FALSE,
                    your_title = "",
                    grid = TRUE,
                    label_size = 10,
                    plot_type = "color",
                    no_negatives = FALSE,
                    return_table = FALSE,
                    color_pal = "Blues",
                    number_size = 3,
                    point_scale = 5,
                    minkowski_p = 2) {


  #extracts assay from your_SE
  if (assay %in% names(SummarizedExperiment::assays(your_SE)) == FALSE){
    stop("The specified assay is not found in your_SE.")
  }

  plotting_data <- SummarizedExperiment::assays(your_SE)[[assay]]

  if(ncol(plotting_data) < 2){
    stop("your_SE must contain at least 2 samples (columns).")
  }

  plotting_data_columns <- colnames(plotting_data)

  if(!(tolower(dist_method) %in% tolower(proxy::pr_DB$get_entry_names()))){
    stop(paste0("\"", dist_method,"\" is not a valid distance/similarity measure. Please choose a similarity or distance measure from those specified in the help page or in summary(proxy::pr_DB)"))
  }


  is_distance_measure <- proxy::pr_DB$get_entry(name = dist_method)[["distance"]]
  if(is_distance_measure){
    if(tolower(dist_method) == "minkowski"){
      pairwise_object <- proxy::dist(plotting_data, method = dist_method, diag = TRUE, by_rows = FALSE, p = minkowski_p)
    } else{
      pairwise_object <- proxy::dist(plotting_data, method = dist_method, diag = TRUE, by_rows = FALSE)
    }
    pairwise_data <- as.matrix(pairwise_object)
    diag(pairwise_data) <- 0
    measure_string <- paste0(dist_method, " distance")
    color_scale <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, color_pal)))(255) #color scale snippet from DESeq2 vignette
    if(cluster_tree){
      sample_hclust <- hclust(pairwise_object)
      sample_levels <- sample_hclust$labels[sample_hclust$order]
      plot_labels <- plot_labels[match(sample_levels, plotting_data_columns)]
    } else {
      sample_levels <- plotting_data_columns
    }
  } else { #is similarity measure
    pairwise_object <- proxy::simil(plotting_data, method = dist_method, diag = TRUE, by_rows = FALSE)
    pairwise_data <- as.matrix(pairwise_object)
    diag(pairwise_data) <- 1
    measure_string <- paste0(dist_method, " similarity")
    color_scale <- colorRampPalette(RColorBrewer::brewer.pal(9, color_pal))(255) #color scale snippet from DESeq2 vignette
    if(cluster_tree){
      pairwise_dist <- proxy::as.dist(pairwise_object)
      sample_hclust <- hclust(pairwise_dist)
      sample_levels <- sample_hclust$labels[sample_hclust$order]
      plot_labels <- plot_labels[match(sample_levels, plotting_data_columns)]
    } else {
      sample_levels <- plotting_data_columns
    }
  }



  pairwise_data <- as.data.frame(pairwise_data)
  pairwise_data_long <- pairwise_data %>%
    tibble::rownames_to_column(var = "sample_i") %>%
    tidyr::pivot_longer(-sample_i, names_to = "sample_j", values_to = "measure") %>%
    dplyr::mutate(sample_i = factor(sample_i, levels = sample_levels)) %>%
    dplyr::mutate(sample_j = factor(sample_j, levels = rev(sample_levels)))

  if(return_table){
    return(pairwise_data_long)
  }

  gg_distplot <- ggplot2::ggplot(pairwise_data_long, ggplot2::aes(x = sample_i, y = sample_j)) +
    ggplot2::scale_x_discrete(position = "top", labels = plot_labels) +
    ggplot2::scale_y_discrete(labels = rev(plot_labels)) +
    ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                   rect = ggplot2::element_blank(),
                   text = ggplot2::element_text(size = label_size),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5),
                   axis.title = ggplot2::element_blank())+
    ggplot2::ggtitle(your_title)

  if(plot_type == "color"){
    gg_distplot <- gg_distplot +
      ggplot2::geom_tile(ggplot2::aes(fill = measure), color = "black") +
      ggplot2::scale_fill_gradientn(colors = color_scale, name = measure_string)
  } else if (plot_type == "circle"){
    gg_distplot <- gg_distplot +
      ggplot2::geom_tile(color = ifelse(grid, "black", "white"), fill = "white") +
      ggplot2::geom_point(ggplot2::aes(size = abs(measure), fill = measure), shape = 21)+
      ggplot2::scale_size(paste0("|", measure_string,"|"), range = c(0, point_scale)) +
      ggplot2::scale_fill_gradientn(colors = color_scale, name = measure_string)
  } else if (plot_type == "number") {
    gg_distplot <- gg_distplot +
      ggplot2::geom_tile(ggplot2::aes(fill = measure), color = "black") +
      ggplot2::scale_fill_gradientn(colors = color_scale, name = measure_string)+
      ggplot2::geom_text(ggplot2::aes(label = round(measure, digits = 2)), color = "black", size = number_size)
  } else {
    stop("plot_type must be one of \"color\", \"circle\", or \"number\"")
  }

  if(cluster_tree){
    dendro_data <- ggdendro::dendro_data(sample_hclust, type = 'rectangle')
    gg_dendrogram <- ggplot2::ggplot(ggdendro::segment(dendro_data))+
      ggplot2::geom_segment(ggplot2::aes(x = x, y = y, xend = xend, yend = yend))+
      ggplot2::ylab(NULL)+
      ggplot2::xlab(NULL)+
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = label_size),
        axis.text = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = "white",colour = "white"),
        axis.ticks = ggplot2::element_blank()
      )
    left_gg_dendogram <- gg_dendrogram + ggplot2::coord_flip()+
      ggplot2::scale_x_reverse(expand = c(0,0.5))+
      ggplot2::scale_y_reverse(expand = c(0,0))
    gg_distplot <- cowplot::plot_grid(left_gg_dendogram, gg_distplot, rel_widths = c(1,7), align = "h", axis = "bt")

  }

  return(gg_distplot)

}
