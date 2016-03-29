#'@title BCheatmap (Barcode Heatmap)
#'
#'@description Creates a heatmap using the top 'n' rows from each column. More or less a wrapper for gplots' "heatmap.2" function :)
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param names Names for each column of the data frame.
#'@param n_clones The top 'n' clones to plot.
#'@param your_title The title for the plot.
#'@param grid Logical. Include a grid or not in the heatmap.
#'@param columnLabels The size of the column labels.
#'@param dendro Logical. Whether or not to show row dendrogram for hierarchical clustering.
#'@param star_size The numerical size of the cell note labels.
#'@param printtable Logical. Prints percent contribution as a table instead of plotting it.
#'@param table_option Character. One of "logs", "reads", or "percents".
#'@param log_transform Logical. Log transform data before clustering and plotting.
#'@param log_choice Data is log transformed with this log.
#'@param variable_log_min Logical. TRUE to make log(0) = minimum nonzero - 1. FALSE to make log(0) = log(100/4000000) - 1.
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

BCheatmap <- function(your_data, names = colnames(your_data), n_clones = 10,
                      your_title = "", grid = TRUE, columnLabels = 1, dendro = FALSE,
                      star_size = 1, printtable = FALSE,
                      table_option = "percents", log_transform = TRUE, log_choice = exp(1),
                      variable_log_min = TRUE,
                      distance_method = "Euclidean", minkowski_power = 1,
                      cellnote_option = "stars", hclust_linkage = "complete",
                      row_order = "hierarchical", clusters = 3) {

  #scales all data to be a percentage of reads instead of number of reads and keeps copy of raw read number
  your_data_list <- list(raw_reads = your_data, prop_table = as.data.frame(prop.table(as.matrix(your_data),2)))
  your_data_list$prop_table[is.na(your_data)] <- 0
  if (any(colSums(your_data_list$prop_table) != 1)){
    stop("One of your columns contained no data")
  }
  your_data_list$prop_table[your_data_list$prop_table == 0] <- NA

  #creates data frame that shows rank of original your_data
  your_data_list <- c(your_data_list, list(ranks = apply(-your_data_list$prop_table, 2, rank, ties.method = "min", na.last = "keep")))

  #subsets those barcodes that have at least one top N clone
  top_clones_choices <- apply(your_data_list$ranks, 1, function(x){any(x<=n_clones, na.rm = TRUE)})
  your_data_list <- lapply(your_data_list, function(x) {x[top_clones_choices,]})

  #creates data frame with '*' for those cells w/ top clones
  your_data_list <- c(your_data_list, list(cellnote_matrix = your_data_list$ranks))
  your_data_list$cellnote_matrix[your_data_list$cellnote_matrix > n_clones] <- NA
  your_data_list$cellnote_matrix[your_data_list$cellnote_matrix <= n_clones] <- "*"
  your_data_list$prop_table[is.na(your_data_list$prop_table)] <- 0

  #takes log of data
  your_data_list <- c(your_data_list, list(logged_data = custom_log(your_data_list$prop_table, log_choice, variable_log_min)))

  #This does the clustering



  if(row_order == "hierarchical") {
    if(distance_method == "Minkowski"){
      hclustering <-hclust(proxy::dist((if (log_transform) your_data_list$logged_data else your_data_list$prop_table), method = distance_method, p = minkowski_power), method = hclust_linkage)
    } else {
      hclustering <-hclust(proxy::dist((if (log_transform) your_data_list$logged_data else your_data_list$prop_table), method = distance_method), method = hclust_linkage)
    }
    e <- hclustering$order
    if(clusters > 0){
      myPalette <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
                     "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
                     "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D",
                     "#8A7C64", "#599861")
      cuts <- cutree(hclustering, k = clusters)
      cluster_colors = myPalette[1:clusters][cuts]
    }
  } else if(row_order == "emergence") {
    e <- do.call(order, -as.data.frame(your_data_list$prop_table))
  } else {
    stop("row_order must be one of \" hierarchical\" or \"emergence\"")
  }



  if(printtable == TRUE){
    switch(table_option,
           reads = return(your_data_list$raw_reads[e,]),
           percents = return(your_data_list$prop_table[e,]),
           logs = return(your_data_list$logged_data[e,]),
           ranks = return(your_data_list$ranks[e,]))
  } else {

    if (row_order == "emergence") {
      call_list <- list(
        x = as.matrix(if(log_transform) your_data_list$logged_data[e,] else your_data_list$prop_table[e,]),
        cellnote = switch(cellnote_option,
                          stars = your_data_list$cellnote_matrix[e,],
                          reads = your_data_list$raw_reads[e,],
                          percents = t(apply(round(your_data_list$prop_table *100, digits = 2), 1, paste0, "%"))[e,],
                          logs = round(your_data_list$logged_data, digits = 2)[e,],
                          ranks = your_data_list$ranks[e,]),
        dendrogram = "none",
        Rowv = FALSE)
    } else if (row_order == "hierarchical") {
      call_list <- list(
        x = as.matrix(if(log_transform) your_data_list$logged_data else your_data_list$prop_table),
        cellnote = switch(cellnote_option,
                          stars = your_data_list$cellnote_matrix,
                          reads = your_data_list$raw_reads,
                          percents = t(apply(round(your_data_list$prop_table *100, digits = 2), 1, paste0, "%")),
                          logs = round(your_data_list$logged_data, digits = 2),
                          ranks = your_data_list$ranks),
        Rowv = as.dendrogram(hclustering),
        dendrogram = if(dendro) "row" else "none"
      )
      if (clusters > 0) {
        call_list <- c(call_list, list(RowSideColors = cluster_colors))
      }
    } else {
      stop("row_order must be one of \" hierarchical\" or \"emergence\"")
    }


    call_list <- c(call_list, list(
      labRow = "",
      notecex = star_size,
      sepwidth=c(0,0.0001),
      sepcolor="black",
      colsep= if (grid) 1:ncol(your_data_list$raw_reads) else NULL,
      rowsep= if (grid) 1:nrow(your_data_list$raw_reads) else NULL,
      offsetRow = -0.4,
      offsetCol = 0,
      notecol = "black",
      margins = c(15,6),
      density.info="none",
      trace = "none",
      symm=F,
      scale="none",
      Colv = names,
      col=(rainbow(256, s = 1, v = 1, start = 0, end = 0.75, alpha = 1)),
      cexCol = columnLabels,
      labCol = names,
      symkey = FALSE,
      key = TRUE,
      keysize = 0.8,
      srtCol = 45,
      main = paste0("\n\n", your_title)
    ))
    do.call(gplots::heatmap.2, call_list)


  }
}

custom_log <- function(x, log_choice, vary_log_min){
  x <- log(x, log_choice)
  #sets the -inf values to be the minimum of the possible value in the data frame -1
  if(vary_log_min) {
    x[x == -Inf] <- (min(x[x > -Inf]) - 1)
  } else {
    x[x == -Inf] <- (log(100/4000000) - 1)
  }
  return(x)
}

