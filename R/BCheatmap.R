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
#'@param dendro Whether or not to rearrange data with dendrogram.
#'@param star_size The numerical size of the cell note labels.
#'@param printtable Logical. Prints percent contribution as a table instead of plotting it.
#'@param table_option Character. One of "logs", "reads", or "percents".
#'@param log_transform Logical. Log transform data before clustering and plotting.
#'@param log_choice Data is log transformed with this log.
#'@param distance_method Character. One of "euclidean", "minkowski", "canberra", manhatan", "maximum", or "binary"
#'@param minkowski_power The power of the Minkowski distance (if minkowski is used).
#'@param cellnote_option Character. One of "stars", "reads", or "percents"
#'@param hclust_linkage Character. One of one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#'@return Displays a heatmap in the current plot window.
#'@examples
#'BCheatmap(your_data = zh33, names = colnames(zh33), n_clones = 10,
#'       your_title = "First Time Point", grid = TRUE, columnLabels = 3)
#'BCheatmap(your_data = zh33, n_clones = 10, printtable = TRUE)
#'@export


BCheatmap2 <- function(your_data, names = colnames(your_data), n_clones = 10,
                      your_title = "", grid = TRUE, columnLabels = 1, dendro = "none",
                      star_size = 1, printtable = FALSE,
                      table_option = "percents", log_transform = TRUE, log_choice = exp(1),
                      distance_method = "euclidean", minkowski_power = 1,
                      cellnote_option = "stars", hclust_linkage = "complete") {

  #scales all data to be a percentage of reads instead of number of reads and keeps copy of raw read number
  raw_reads <- your_data
  your_data <- as.data.frame(prop.table(as.matrix(your_data),2))
  your_data[is.na(your_data)] <- 0
  if (any(colSums(your_data) != 1)){
    stop("One of your columns contained no data")
  }
  your_data[your_data == 0] <- NA

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

  your_data[is.na(your_data)] <- 0

  #takes log of data
  your_data_log <- custom_log(your_data, log_choice)
  e = hclust(dist((if (log_transform) your_data_log else your_data), method = distance_method, p = minkowski_power),
             method = hclust_linkage)$order

  if(printtable == TRUE){
    switch(table_option,
           reads = return(raw_reads[e,]),
           percents = return(your_data[e,]),
           logs = return(your_data_log[e,]))
  } else {
    gplots::heatmap.2(as.matrix((if (log_transform) your_data_log[e,] else your_data[e,])),
                      labRow = "",
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

                      key = TRUE,
                      keysize = 0.8,
                      srtCol = 45,
                      main = your_title)
  }
}

custom_log <- function(x, log_choice, inf_scaler){
  x <- log(x, log_choice)
  #sets the -inf values to be the minimum of the possible value in the data frame -1
  x[x == -Inf] <- (min(x[x > -Inf]) - 1)
  return(x)
}

