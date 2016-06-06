#'@title Barcode Log Read Fraction of Hematopoiesis
#'
#'@description Creates a heatmap using the top 'n' rows from each column. Data undergoes a natural log transformation before being plotted. Optionally, will write a percentage key to the current directory.
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param names Names for each column of the data frame.
#'@param n_clones The top 'n' clones to plot.
#'@param your_title The title for the plot.
#'@param grid Logical. Include a grid or not in the heatmap.
#'@param columnLabels The size of the column labels.
#'@param dendro Whether or not to rearrange data with dendrogram.
#'@param printtable Logical. Prints percent contribution as a table instead of plotting it.
#'@param log_transorm Logical. Log transform the data or not.
#'@param log_choice Data is log transformed with this log
#'@return Displays a heatmap in the current plot window.
#'@examples
#'BCLH(your_data = zh33, names = colnames(zh33), n_clones = 10,
#'       your_title = "First Time Point", grid = TRUE, columnLabels = 3)
#'BCLH(your_data = zh33, n_clones = 10, printtable = TRUE)
#'@export


BCLH <- function(your_data, names = colnames(your_data), n_clones = 10, your_title = "", grid = TRUE, columnLabels = 1, dendro = "none", star_size = 1, printtable = FALSE, log_choice = exp(1)) {

  #scales all data to be a percentage of reads instead of number of reads
  your_data <- as.data.frame(prop.table(as.matrix(your_data),2))
  your_data[is.na(your_data)] <- 0


  #this gets the top clones from each of the samples, and then the following 'unique()' function
  #will combine this list
  order.list <- apply(your_data, 2, function(x){sort.list(-x)[1:n_clones]})
  top_clones = unique(as.vector(order.list))

  #takes only the top clones subset of the total data
  your_data <- your_data[top_clones,]

  #creates a copy of the data
  list_topclones <- your_data

  #orders each column according to the top clones of each
  for(i in 1:ncol(your_data)){
    list_topclones[,i] <- order(-your_data[,i])
  }

  #creates empty data frame with dimension of subsetted data
  is_a_topclone <- data.frame(matrix(ncol = ncol(your_data), nrow = nrow(your_data)))

  #populates the empty data frame with '*' for those cells w/ top clones
  for(i in 1:ncol(your_data)){
    is_a_topclone[list_topclones[1:n_clones,i],i] <- "*"
  }

  #takes the log of the data set
  data.log <- log(your_data, log_choice)

  #sets the -inf values to be the minimum of the possible value in the data frame -1
  data.log[data.log == -Inf] <- (min(data.log[data.log > -Inf]) - 1)

  e = hclust(dist(data.log))$order

  if(grid == TRUE){
    columnsep <- 1:ncol(data.log)
    rowseparation <- 1:nrow(data.log)
  } else {
    columnsep <- NULL
    rowseparation <- NULL
  }

  if(printtable == TRUE){
    return(your_data[e,])
  } else {
    gplots::heatmap.2(as.matrix(data.log[e,]),
                      labRow = "",
                      notecex = star_size,
                      sepwidth=c(0,0.0001),
                      sepcolor="black",
                      colsep=columnsep,
                      rowsep=rowseparation,
                      offsetRow = -0.4,
                      offsetCol = 0,
                      Rowv = FALSE,
                      cellnote = is_a_topclone[e,], #applies the stars to the correct cells according to the hclust 'e'
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
