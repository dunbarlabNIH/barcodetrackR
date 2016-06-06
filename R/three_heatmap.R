three_heatmap <- function(your_data, names = colnames(your_data), n_clones = 10,
                          your_title = "", grid = TRUE, columnLabels = 1, dendro = "none",
                          star_size = 1, printtable = FALSE,
                          table_option = "percents", log_transform = TRUE, log_choice = exp(1),
                          variable_log_min = FALSE,
                          distance_method = "Euclidean", minkowski_power = 1,
                          cellnote_option = "stars", hclust_linkage = "complete",
                          row_order = "hierarchical", clusters = 3) {
  #scales all data to be a percentage of reads instead of number of reads and keeps copy of raw read number
  raw_reads <- your_data
  your_data <- as.data.frame(prop.table(as.matrix(your_data),2))
  your_data[is.na(your_data)] <- 0
  if (any(colSums(your_data) != 1)){
    stop("One of your columns contained no data")
  }
  your_data[your_data == 0] <- NA
  lib_IDs <- unique(substring(rownames(your_data),1,6))
  data_list <- list()
  plot_list <- list()
  for(i in 1:length(lib_IDs)) {
        data_list[[i]] <- your_data[substring(rownames(your_data),1,6) == lib_IDs[i],]
  }
  
  grab_grob <- function(){
    gridGraphics::grid.echo()
    grid::grid.grab()
  }
  
  gl <- lapply(1:3, function(i){
    small_heatmap(data_list[[i]], names, n_clones, your_title = lib_IDs[[i]],
                  grid,
                  columnLabels,
                  dendro,
                  star_size,
                  printtable,
                  table_option, log_transform, log_choice,
                  variable_log_min,
                  distance_method, minkowski_power,
                  cellnote_option, hclust_linkage,
                  row_order, clusters)
    grab_grob()
  })
  
  
  grid::grid.newpage()
  gridExtra::grid.arrange(grobs=gl, ncol=2, clip=TRUE)
  
}

custom_log_II <- function(x, log_choice, vary_log_min){
  x <- log(x, log_choice)
  #sets the -inf values to be the minimum of the possible value in the data frame -1
  if(vary_log_min) {
    x[x == -Inf] <- (min(x[x > -Inf]) - 1)
  } else {
    x[x == -Inf] <- (log(100/4000000) - 1)
  }
  return(x)
}




small_heatmap <- function(your_data, names = colnames(your_data), n_clones = 10,
                          your_title = "", grid = TRUE, columnLabels = 1, dendro = "none",
                          star_size = 1, printtable = FALSE,
                          table_option = "percents", log_transform = TRUE, log_choice = exp(1),
                          variable_log_min = TRUE,
                          distance_method = "Euclidean", minkowski_power = 1,
                          cellnote_option = "stars", hclust_linkage = "complete",
                          row_order = "hierarchical", clusters = 3) {
  
  raw_reads <- your_data
  #creates data frame that shows rank of original your_data
  your_data_ranked <- apply(-your_data, 2, rank, ties.method = "min", na.last = "keep")
  
  #subsets those barcodes that have at least one top N clone
  top_clones_choices <- apply(your_data_ranked, 1, function(x){any(x<=n_clones, na.rm = TRUE)}) 
  your_data <- your_data[top_clones_choices,]
  raw_reads <- raw_reads[top_clones_choices,]
  
  
  
  #creates empty data frame with dimension of subsetted data &
  #populates the empty data frame with '*' for those cells w/ top clones
  your_data_ranked <- apply(-your_data, 2, rank, ties.method = "min", na.last = "keep")
  is_a_topclone <- data.frame(matrix(ncol = ncol(your_data), nrow = nrow(your_data)))
  is_a_topclone[your_data_ranked <= n_clones] <- "*"
  is_a_topclone[your_data_ranked > n_clones] <- ""
  is_a_topclone[is.na(is_a_topclone)] <- ""
  
  your_data[is.na(your_data)] <- 0
  
  #takes log of data
  your_data_log <- custom_log_II(your_data, log_choice, variable_log_min)

  
  if(row_order == "hierarchical") {
    
    if(distance_method == "Minkowski"){
      hclustering <-hclust(proxy::dist((if (log_transform) your_data_log else your_data), method = distance_method, p = minkowski_power), method = hclust_linkage)
    } else {
      hclustering <-hclust(proxy::dist((if (log_transform) your_data_log else your_data), method = distance_method), method = hclust_linkage)
    }
    e <- hclustering$order
    if( clusters > 0){
      cuts <- cutree(hclustering, k = clusters)
      cluster_colors = RColorBrewer::brewer.pal(clusters, "Set1")[cuts[e]]
    } else {
      cluster_colors = rep("white", length(e))
    }
    
  } else if(row_order == "emergence") {
    e <- do.call(order, -as.data.frame(your_data_log))
  } else {
    stop("row_order must be one of \" hierarchical\" or \"emergence\"")
  }
  
  if(printtable == TRUE){
    switch(table_option,
           reads = return(raw_reads[e,]),
           percents = return(your_data[e,]),
           logs = return(your_data_log[e,]))
  } else {
    gplots::heatmap.2(as.matrix((if (log_transform) your_data_log[e,] else your_data[e,])),
                      labRow = rep("", nrow(your_data_log)),
                      notecex = star_size,
                      sepwidth=c(0,0.0001),
                      sepcolor="black",
                      colsep= if (grid) 1:ncol(your_data) else NULL,
                      rowsep= if (grid) 1:nrow(your_data) else NULL,
                      offsetRow = -0.4,
                      offsetCol = 0,
                      Rowv = FALSE,
                      cellnote = switch(cellnote_option,
                                        stars = is_a_topclone[e,],
                                        reads = (raw_reads[rownames(your_data[e,]),]),
                                        percents = t(apply(round(your_data[e,]*100, digits = 2), 1, paste0, "%")),
                                        logs = round(your_data_log[e,], digits = 2)),
                      dendrogram = dendro,
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
                      RowSideColors = if(row_order == "hierarchical") cluster_colors else rep("white", length(e)),
                      key = TRUE,
                      keysize = 0.8,
                      srtCol = 45,
                      main =your_title)
  }
}


