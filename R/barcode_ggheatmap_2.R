#'@importFrom proxy dist
#'@importFrom SummarizedExperiment assays
#'
#'@title barcode_ggheatmap (Barcode Heatmap using ggplot2)
#'
#'@description Creates a heatmap using the top 'n' rows from each column, using ggplot2!
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param hclust_assay Perform hierarchical clustering on this assay.
#'@param visual_assay Use this assay for visualization.
#'@param names Vector of x axis labels.
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
#'@param row_order Character. One of "hierarchical" or "emergence" to organize rows.
#'@param clusters How many clusters to cut hierarchical tree into for display when row_order is "hierarchical".
#'@return Displays a heatmap in the current plot window.
#'@examples
#'BCheatmap(your_data = zh33, names = colnames(zh33), n_clones = 10,
#'       your_title = "First Time Point", grid = TRUE, columnLabels = 3)
#'BCheatmap(your_data = zh33, n_clones = 10, printtable = TRUE)
#'@export
barcode_ggheatmap_2 <- function(your_SE,
                                hclust_assay = "logs",
                                visual_assay = "logs",
                                selections = NULL,
                                names = colData(your_SE),
                                n_clones = 10,
                                your_title = "",
                                grid = TRUE,
                                label_size = 1,
                                dendro = FALSE,
                                cellnote_size = 4,
                                printtable = FALSE,
                                table_option = "percents",
                                log_transform = TRUE,
                                distance_method = "Euclidean",
                                minkowski_power = 1,
                                cellnote_option = "stars",
                                hclust_linkage = "complete",
                                row_order = "hierarchical",
                                clusters = 0) {

  #subset your SE based on your selections in selections (if any)
  your_SE <- subset_SE(your_SE, selections)
  #subsets those barcodes that have at least one top N clone
  top_clones_choices <- apply(assays(your_SE)$ranks, 1, function(x){any(x<=n_clones, na.rm = TRUE)})
  your_SE <- your_SE[top_clones_choices,]

  #creates data frame with '*' for those cells w/ top clones
  cellnote_matrix = assays(your_SE)$ranks
  cellnote_matrix[cellnote_matrix > n_clones] <- NA
  cellnote_matrix[cellnote_matrix <= n_clones] <- "*"

  if(row_order == "hierarchical") {
    if(distance_method == "Minkowski"){
      hclustering <- hclust(proxy::dist((if (log_transform) assays(your_SE)$logs else assays(your_SE)$percentages), method = distance_method, p = minkowski_power), method = hclust_linkage)
    } else {
      hclustering <- hclust(proxy::dist((if (log_transform) assays(your_SE)$logs else assays(your_SE)$percentages), method = distance_method), method = hclust_linkage)
    }
    barcode_order <- rev(hclustering$order)
  } else if(row_order == "emergence") {
    barcode_order <- do.call(order, assays(your_SE)$percentages)
  } else {
    stop("row_order must be one of \" hierarchical\" or \"emergence\"")
  }




  if(v){
    plotting_data <- reshape2::melt(as.matrix(your_data_list$logged_data))
    actual_scale <- c(log(100/4000000, log_choice) - 1, log(100/4000000, log_choice), log(0.001, log_choice), log(0.01,log_choice), log(0.1, log_choice), 0)
    your_scale <- scales::rescale(actual_scale, to = c(0,1))
    your_breaks <- c(actual_scale)
    your_labels <- c("Defined 0", "0.0025% ", "0.1%", "1%", "10%", "100%")
    your_limits <- c(log(100/4000000, log_choice)-1,0)
    your_colors <- c("#4575B4", "#4575B4", "lightblue", "#fefeb9", "#D73027", "red4")
  } else {
    plotting_data <- reshape2::melt(as.matrix(your_data_list$prop_table))
    your_scale <- c(0, 100/4000000, 0.005, 0.01, 0.1, 1)
    your_breaks <- c(0, 100/4000000, 0.005, 0.01, 0.1, 1)
    your_labels <- c("0%", "0.0025%", "0.5%", "1%", "10%", "100%")
    your_limits <- c(0,1)
    your_colors <- c("#4575B4", "#4575B4", "lightblue", "#fefeb9", "#D73027", "red4")
  }



  colnames(plotting_data) <- c("BARCODE", "SAMPLE", "Size")
  plotting_data$BARCODE <- factor(plotting_data$BARCODE, levels = rev(rownames(your_data_list$prop_table)[barcode_order]))

  if (cellnote_option %in% c("logs", "ranks", "stars", "percents", "reads")){
    if(cellnote_option == "logs"){
      plotting_data$CELLNOTE <- round(reshape2::melt(as.matrix(your_data_list$logged_data))$value, 2)
      cellnote_vjust = 0.5
    }
    if(cellnote_option == "ranks"){
      plotting_data$CELLNOTE <- reshape2::melt(as.matrix(your_data_list$ranks))$value
      cellnote_vjust = 0.5
    }
    if(cellnote_option == "stars"){
      plotting_data$CELLNOTE <- reshape2::melt(as.matrix(your_data_list$cellnote_matrix))$value
      cellnote_vjust = 0.75
    }
    if(cellnote_option == "percents"){
      plotting_data$CELLNOTE <- paste0(round(reshape2::melt(as.matrix(your_data_list$prop_table))$value*100,2), "%")
      cellnote_vjust = 0.5
    }
    if(cellnote_option == "reads"){
      plotting_data$CELLNOTE <- reshape2::melt(as.matrix(your_data_list$raw_reads))$value
      cellnote_vjust = 0.5
    }
  } else {
    stop("cellnote_option must be one of c(\"logs\", \"ranks\", \"stars\", \"percents\", \"reads\")")
  }

  if(grid)gridColor = "black" else gridColor = NA


  false_column_label <- unique(colnames(your_data_list$logged_data))[which(max(nchar(unique(colnames(your_data_list$logged_data)))) == nchar(unique(colnames(your_data_list$logged_data))))[1]]

  if(printtable == TRUE){
    switch(table_option,
           reads = return(your_data_list$raw_reads[barcode_order,]),
           percents = return(your_data_list$prop_table[barcode_order,]),
           logs = return(your_data_list$logged_data[barcode_order,]),
           ranks = return(your_data_list$ranks[barcode_order,]))
  } else {

    g1_heatmap <- ggplot2::ggplot(plotting_data, ggplot2::aes(x = SAMPLE, y = BARCODE))+
      ggplot2::geom_tile(ggplot2::aes(fill = Size), color = gridColor)+
      ggplot2::geom_text(ggplot2::aes(label = CELLNOTE), vjust = cellnote_vjust, size = cellnote_size)+
      ggplot2::scale_fill_gradientn(
        colors = your_colors,
        values = your_scale,
        breaks = your_breaks,
        labels = your_labels,
        limits = your_limits,
        expand = c(0,0))+
      ggplot2::scale_y_discrete(labels = NULL, breaks = NULL, expand = c(0,0))+
      ggplot2::scale_x_discrete(expand = c(0,0), labels = names)+
      ggplot2::ylab(NULL)+
      ggplot2::xlab(NULL)+
      ggplot2::ggtitle(paste0("\n", your_title, "\n"))+
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 20),
        axis.text.x = ggplot2::element_text(angle=90, hjust = 1, vjust = 0.5, size = label_size),
        legend.text = ggplot2::element_text(size =  15, face = 'bold'),
        legend.title = ggplot2::element_text(size =  15),
        axis.ticks = ggplot2::element_blank(),
        legend.key.height=ggplot2::unit(5,"line"))


    if(row_order == 'emergence'){
      g1_heatmap
    } else {
      dendro_data <- ggdendro::dendro_data(hclustering, type = 'rectangle')
      g2_dendrogram <- ggplot2::ggplot(ggdendro::segment(dendro_data))+
        ggplot2::geom_segment(ggplot2::aes(x = x, y = y, xend = xend, yend = yend))+
        ggplot2::scale_x_discrete(expand = c(.5/nrow(your_data_list$logged_data),0.01))+
        ggplot2::scale_y_reverse(expand = c(0.01,0), labels = false_column_label, breaks = if(row_order == 'emergence') NULL else mean(hclustering$height))+
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

      if(dendro){

        if(clusters > 0){
          # custom_colors <- c(RColorBrewer::brewer.pal(8, 'Set1'), "cyan", "black", "grey")
          custom_colors <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588")
          clusters_df <-data.frame(CLUSTERS = cutree(hclustering, clusters), ASSIGNMENT=factor(hclustering$labels,levels=hclustering$labels[(hclustering$order)]))
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
          gridExtra::grid.arrange(g2_dendrogram,
                                  g3_clusters,
                                  g1_heatmap,
                                  ncol = 3,
                                  widths = c(1,.2,4))

        } else {
          gridExtra::grid.arrange(g2_dendrogram,
                                  g1_heatmap,
                                  ncol = 2,
                                  widths = c(1,4))
        }


      } else {
        g1_heatmap
      }
    }




  }

}

