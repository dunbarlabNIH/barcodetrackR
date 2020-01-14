#'@importFrom SummarizedExperiment assays
#'@importFrom rlang %||%
#'@title barcode_ggheatmap (Barcode Heatmap using ggplot2)
#'
#'@description Creates a heatmap using the top 'n' rows from each column in the Summarized Experiment object, using ggplot2.
#'
#'@param your_SE A Summarized Experiment object.
#'@param assay Perform visualization and clustering on this assay.
#'@param plot_labels Vector of x axis labels.
#'@param n_clones The top 'n' clones to plot.
#'@param your_title The title for the plot.
#'@param grid Logical. Include a grid or not in the heatmap.
#'@param label_size The size of the column labels.
#'@param dendro Logical. Whether or not to show row dendrogram for hierarchical clustering.
#'@param cellnote_size The numerical size of the cell note labels.
#'@param printtable Logical. Prints percent contribution as a table instead of plotting it.
#'@param table_option Character. One of "logs", "reads", or "percents" for printing.
#'@param distance_method Character. Use summary(proxy::pr_DB) to see all options.
#'@param minkowski_power The power of the Minkowski distance (if minkowski is used).
#'@param cellnote_option Character. One of "stars", "reads", or "percents"
#'@param hclust_linkage Character. One of one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#'@param row_order Character; "hierarchical" to perform hierarchical clustering on the output and order in that manner, "emergence" to organize rows  by order of presence in data (from left to right), or a character vector of rows within the summarized experiment to plot.
#'@param clusters How many clusters to cut hierarchical tree into for display when row_order is "hierarchical".
#'@return Displays a heatmap in the current plot window.
#'@examples
#'BCheatmap(your_data = zh33, names = colnames(zh33), n_clones = 10,
#'       your_title = "First Time Point", grid = TRUE, columnLabels = 3)
#'BCheatmap(your_data = zh33, n_clones = 10, printtable = TRUE)
#'@export
barcode_ggheatmap_2 <- function(your_SE,
                                assay = "logs",
                                selections = list(),
                                plot_labels = NULL,
                                n_clones = 10,
                                your_title = "",
                                grid = TRUE,
                                label_size = 1,
                                dendro = FALSE,
                                cellnote_size = 4,
                                printtable = FALSE,
                                table_option = "percents",
                                distance_method = "Euclidean",
                                minkowski_power = 2,
                                cellnote_assay = "stars",
                                hclust_linkage = "complete",
                                row_order = "hierarchical",
                                clusters = 0,
                                percent_scale = c(0, 0.000025, 0.001, 0.01, 0.1, 1)) {

  #subset your SE based on your selections in selections (if any)
  if(length(selections) > 0){
    your_SE <- subset_SE(your_SE, selections)
  }

  #get labels for heatmap
  plot_labels <- plot_labels %||% colnames(your_SE)
  if(length(plot_labels) != ncol(your_SE)){
    stop("plot_labels must be same length as number of columns being plotted")
  }


  #subset the rows of the summarized experiment and get the ordering of barcodes within the hematmap for plotting
  if(row_order == "hierarchical" | row_order == "emergence") {

    #subsets those barcodes that have at least one top N clone
    top_clones_choices <- apply(assays(your_SE)$ranks, 1, function(x){any(x<=n_clones, na.rm = TRUE)})
    your_SE <- your_SE[top_clones_choices,]

    #creates data frame with '*' for those cells w/ top clones
    cellnote_matrix = assays(your_SE)$ranks
    cellnote_matrix[cellnote_matrix > n_clones] <- NA
    cellnote_matrix[cellnote_matrix <= n_clones] <- "*"
    assays(your_SE)$stars <- as.data.frame(cellnote_matrix)

    #this does the heavy duty plotting set-up. It sets the order of the data on the heatmap and the dendrogram/cluster cuts
    if(row_order == "hierarchical"){
      clustering_data <- assays(your_SE)[[assay]]
      clustering_data.dist <- proxy::dist(clustering_data, method = distance_method, p = minkowski_power)
      hclustering <- hclust(clustering_data.dist, method = hclust_linkage)
      barcode_order <- rownames(your_SE)[rev(hclustering$order)]

      if(dendro){
        dendro_data <- ggdendro::dendro_data(hclustering, type = 'rectangle')
      }

      if(clusters>0){
        clustercuts_data <-data.frame(clusters = cutree(hclustering, clusters),
                                      assignment = factor(hclustering$labels,levels = hclustering$labels[(hclustering$order)]))
      }


    } else if(row_order == "emergence"){
      barcode_order <- rownames(your_SE)[do.call(order, assays(your_SE)$percentages)]
    }

  } else {
    message("using supplied row_order")
    your_SE <- your_SE[row_order,]
    barcode_order <- row_order
  }

  #create scale for plotting
  log_used <- S4Vectors::metadata(your_SE)$log_base
  scale_factor_used <- S4Vectors::metadata(your_SE)$scale_factor
  percent_scale_labels <- paste0(percent_scale * 100, "%")
  log_scale <- log(percent_scale*scale_factor_used + 1, base = log_used)
  #need to fix this if percent scale is larger than 6 numbers
  color_scale <- c("#4575B4", "#4575B4", "lightblue", "#fefeb9", "#D73027", "red4")

  #organizing data for plotting
  plotting_data <- tibble::rownames_to_column(assays(your_SE)[[assay]], var = "sequence")
  plotting_data <- tidyr::pivot_longer(plotting_data, cols = -sequence, names_to = "sample_name", values_to = "value")
  plotting_data$sample_name <- factor(plotting_data$sample_name, levels = names)
  plotting_data$sequence <- factor(plotting_data$sequence, levels = barcode_order)

  #organizing labels for plotting overlay
  plotting_cellnote <- tibble::rownames_to_column(assays(your_SE)[[cellnote_assay]], var = "sequence")
  plotting_cellnote <- tidyr::pivot_longer(plotting_cellnote, cols = -sequence, names_to = "sample_name", values_to = "label")
  plotting_data$cellnote <- plotting_cellnote$label
  if(is.numeric(plotting_data$cellnote)){
    if(cellnote_assay == "percentages"){
      plotting_data$cellnote <- paste0(round(plotting_data$cellnote*100, digits = 2), "%")
    } else {
      plotting_data$cellnote <- round(plotting_data$cellnote, digits = 2)
    }
  }

  if(grid) grid_color = "black" else grid_color = NA

  g1_heatmap <- ggplot2::ggplot(plotting_data, ggplot2::aes(x = sample_name, y = sequence))+
    ggplot2::geom_tile(ggplot2::aes(fill = value), color = grid_color)+
    ggplot2::geom_text(ggplot2::aes(label = cellnote), vjust = 0.75, size = cellnote_size, color = "black")+
    ggplot2::scale_fill_gradientn(
      colors = color_scale,
      values = scales::rescale(log_scale, to = c(0,1)),
      breaks = log_scale,
      limits = c(min(log_scale), max(log_scale)),
      labels = percent_scale.labels,
      expand = c(0,0)
    )+
    ggplot2::scale_y_discrete(labels = NULL, breaks = NULL, expand = c(0,0))+
    ggplot2::scale_x_discrete(expand = c(0,0))+
    ggplot2::ylab(NULL)+
    ggplot2::xlab(NULL)+
    ggplot2::ggtitle(paste0("\n", your_title, "\n"))+
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 20),
      axis.text.x = ggplot2::element_text(angle=90, hjust = 1, vjust = 0.5, size = label_size),
      #legend.text = ggplot2::element_text(size =  15, face = 'bold'),
      legend.title = ggplot2::element_text(size =  15),
      axis.ticks = ggplot2::element_blank())

  if(row_order == 'emergence'){
    g1_heatmap
  } else {
    if(dendro){
      g2_dendrogram <- ggplot2::ggplot(ggdendro::segment(dendro_data))+
        ggplot2::geom_segment(ggplot2::aes(x = x, y = y, xend = xend, yend = yend))+
        ggplot2::scale_x_discrete(expand = c(.5/nrow(your_SE),0.01))+
        ggplot2::scale_y_reverse(expand = c(0.01,0), labels = "false_column_label", breaks = if(row_order == 'emergence') NULL else mean(hclustering$height))+
        ggplot2::coord_flip()+
        ggplot2::ylab(NULL)+
        ggplot2::xlab(NULL)+
        ggplot2::ggtitle("\n \n")+
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 20),
          axis.text.x = ggplot2::element_text(colour = 'white', angle=90, hjust = 1, vjust = 0.5, size = label_size),
          panel.background = ggplot2::element_rect(fill = "white",colour = "white"),
          axis.ticks = ggplot2::element_blank()
        )
    }
    if(clusters > 0){
      g3_clusters <-ggplot2::ggplot(clusters_df, ggplot2::aes(x =1,y=ASSIGNMENT,fill=factor(CLUSTERS)))+
        ggplot2::geom_tile()+
        ggplot2::ggtitle("\n \n")+
        ggplot2::scale_fill_manual(values = custom_colors[1:clusters])+
        ggplot2::scale_x_continuous(expand=c(0,0), labels = false_column_label, breaks = 1)+
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 20),
          axis.title=ggplot2::element_blank(),
          axis.ticks=ggplot2::element_blank(),
          axis.text.y=ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(colour = 'white', angle=90, hjust = 1, vjust = 0.5, size = label_size),
          legend.position="none")
    }
  }

  #now finally plot using cowplot
  if(clusters > 0 & dendro){
    cowplot::plot_grid(g2_dendrogram, g3_clusters, g1_heatmap, rel_widths = c(1, .2, 4), ncol = 3)
  } else if (clusters == 0 & dendro){
    cowplot::plot_grid(g2_dendrogram, g1_heatmap, rel_widths = c(1, 4), ncol = 2)
  } else if (clusters > 0 & !dendro){
    cowplot::plot_grid(g3_clusters, g1_heatmap, rel_widths = c(1, .2), ncol = 2)
  } else if (clusters == 0 & !dendro){
    g1_heatmap
  }

}

