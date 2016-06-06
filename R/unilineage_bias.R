#'@title unilineage_bias
#'
#'@description Methods for displaying clones biased towards a specific lineage.
#'@param your_data A data frame. Arranged by time point.
#'@param CT The number of cell types in your data frame.
#'@param TP The number of time points in your data frame.
#'@param percent_thresh The percent threshold the unilineage clones must be above.
#'@param ratio_thresh The ratio unilineage clones must be above to be considered biased.
#'@param line_months Numerical vector giving the exact months of each of the time points.
#'@param only_biased Logical. Whether or not to plot only the biased clones.
#'@param plot_mode Character. One of "heatmap", "bars", or "lines".
#'@param text_size Numerical size of all the text.
#'@param line_size Numerical size of the line thickness in "lines".
#'@param dot_size Numerical size of the dots in "lines".
#'@param your_title Character. Your title for the plot.
#'@param y_upper Numerical y upper limit for the "lines" and "bars" plots.
#'@param y_lower Numerical y lower limit for the "lines" and "bars" plots.
#'@param cellnote_display Character. One of "stars" or "percents" for "heatmap" plots.
#'@param by_celltype Logical. Whether or not to order by cell type instead of time point.
#'@param print_table Logical. Prints the table to be plotted rather than plotting it.
#'@return Displays a heatmap, line plot, or bar chart in the current plot window.
#'@examples
#'unilineage_bias(your_data = zh33_byTP, CT = 5, TP = 7, your_title = "NK Bias", plot_mode = "heatmap")
#'@export
unilineage_bias <- function(your_data, CT, TP, percent_thresh = 0.01,
                            ratio_thresh = 10, line_months = c(1:TP),
                            only_biased = FALSE, plot_mode = "",
                            text_size = 20, line_size = 3,
                            dot_size = 4, your_title = "",
                            y_upper = 1, y_lower = 0,
                            cellnote_display = "stars", by_celltype = FALSE,
                            print_table = FALSE, cellnote_size = 3){
  print(colnames((your_data)))
  your_data <- as.data.frame(prop.table(as.matrix(your_data), margin = 2))
  your_data[is.na(your_data)] <- 0
  x <- lapply(1:TP, function(i){
    temp <- your_data[,c((CT*i - (CT-1)):(CT*i))]
    temp <- temp[temp[,CT] > percent_thresh,]
    return(rownames(temp[apply(temp[,1:CT-1], 1, max)*ratio_thresh < temp[,CT],]))
  })
  xx <- unique(unlist(x))
  your_data <- your_data[xx,]
  cellnote_matrix <- matrix(nrow = nrow(your_data), ncol = ncol(your_data))
  rownames(cellnote_matrix) <- xx
  for(i in 1:TP){
    cellnote_matrix[rownames(your_data) %in% x[[i]],CT*i] <- "*"
  }
  by_celltype_order <- unlist(lapply(1:CT, function(i){seq(from = i, to = TP*CT, by = CT)}))
  thirtycolors <- c("#9843A1","#CDAB44","#A14854","#8269DD","#DD4469","#A53F2A","#6F8BDD","#7D7078","#C6D949","#6497BF","#80622E","#CD91DC","#D49980","#D74FDB","#E2442C","#76D7D0","#CBC6D0","#5A8A34","#584E8F","#D77E33","#C3D497","#363550","#D6489C","#324026","#6ADA44","#532424","#6AD88C","#578B72","#772B5D","#D087A9")
  thirtycolors <- sample(thirtycolors)

  if(by_celltype){
    total_order <- colnames(your_data)[by_celltype_order]
  } else {
    total_order <- colnames(your_data)
  }

  if(plot_mode == "line"){
    your_data <- your_data[,seq(from = CT, to = CT*TP, by = CT)]
    cellnote_matrix <- cellnote_matrix[,seq(from = CT, to = CT*TP, by = CT)]
    colnames(your_data) <- line_months
    your_data$EMERGENCE <- apply(your_data, 1, which.max)
    your_data$BARCODE <- rownames(your_data)
    melted_data <- reshape2::melt(your_data, id.vars = c("EMERGENCE", "BARCODE"))
    colnames(melted_data) <- c("EMERGENCE", "BARCODE", "SAMPLE", "PROP")
    melted_data$EMERGENCE <- as.character(melted_data$EMERGENCE)
    melted_data$SAMPLE <- as.numeric(as.character(melted_data$SAMPLE))
    if(print_table){
      return(melted_data)
    }
    ggplot2::ggplot(melted_data, ggplot2::aes(x = SAMPLE, y = PROP))+
      ggplot2::scale_x_continuous(breaks = line_months)+
      ggplot2::scale_fill_manual(values = thirtycolors)+
      ggplot2::scale_color_manual(values = thirtycolors)+
      ggplot2::geom_line(ggplot2::aes(group = BARCODE, color = BARCODE, fill = EMERGENCE), size = line_size)+
      ggplot2::geom_point(color = "black", size = dot_size)+
      ggplot2::facet_wrap(~ EMERGENCE)+
      ggplot2::guides(color = FALSE, fill = FALSE)+
      ggplot2::coord_cartesian(ylim = c(y_lower, y_upper))+
      ggplot2::ggtitle(paste0("\n", your_title, "\n"))+
      ggplot2::theme(text = ggplot2::element_text(size=text_size),
                     panel.grid.major.y = ggplot2::element_line(colour = "grey"),
                     axis.ticks.length = ggplot2::unit(1, "cm"),
                     axis.ticks = ggplot2::element_line(size = 1),
                     panel.background = ggplot2::element_rect(fill = "white", colour = "black"),
                     axis.text.x = ggplot2::element_text(angle = 90))
  }else if(plot_mode == "bar"){
    print(colnames(your_data))
    if(only_biased){
      if(ncol(your_data) == CT){
        stop("Cannot use only biased when looking at only one time point")
      }
      your_data <- your_data[,seq(from = CT, to = CT*TP, by = CT)]
      listy <- list()
      for(i in 1:TP){
        listy[[i]] <- reshape2::melt(as.matrix(your_data[x[[i]],i, drop = FALSE]))
      }
    } else {
      listy <- list()
      for(i in 1:TP){
        listy[[i]] <- reshape2::melt(as.matrix(your_data[x[[i]],c((CT*i - (CT-1)):(CT*i)), drop = FALSE]))
      }
    }
    print(listy)
    plotting_table <- do.call(rbind, listy)
    print(plotting_table)
    colnames(plotting_table) <- c("BARCODE", "SAMPLE", "PROP")
    plotting_table$SAMPLE <- factor(plotting_table$SAMPLE, levels = total_order)
    plotting_table <- plotting_table[order(plotting_table$SAMPLE),]
    plotting_table <- plotting_table[!(plotting_table$SAMPLE == "17hm CD56SP FAKE"),]
    if(print_table){
      return(plotting_table)
    }
    ggplot2::ggplot(plotting_table, ggplot2::aes(x = SAMPLE, y = PROP, fill = BARCODE))+
      ggplot2::scale_fill_manual(values = thirtycolors)+
      ggplot2::geom_bar(ggplot2::aes(y = PROP), stat = "identity",  colour = "black")+
      ggplot2::guides(fill = FALSE)+
      ggplot2::scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.1),
                                  labels = function(x){paste0(x*100, "%")},
                                  expand = c(0,0))+
      ggplot2::coord_cartesian(ylim = c(y_lower,y_upper))+
      ggplot2::ggtitle(paste0("\n", your_title, "\n"))+
      ggplot2::theme(text = ggplot2::element_text(size=text_size),
                     panel.grid.major.y = ggplot2::element_line(colour = "grey"),
                     axis.ticks.length = ggplot2::unit(1, "cm"),
                     axis.ticks = ggplot2::element_line(size = 1),
                     panel.background = ggplot2::element_rect(fill = "white", colour = "black"),
                     axis.text.x = ggplot2::element_text(angle = 90))
  }else if(plot_mode == "heatmap"){
    if(only_biased){
      if(ncol(your_data) == CT){
        stop("Cannot use only biased when looking at only one time point")
      }
      your_data <- your_data[,seq(from = CT, to = CT*TP, by = CT)]
      cellnote_matrix <- cellnote_matrix[,seq(from = CT, to = CT*TP, by = CT)]
    } else {
      if(by_celltype){
        your_data <- your_data[,by_celltype_order]
        cellnote_matrix <- cellnote_matrix[,by_celltype_order]
      }
    }

    if(cellnote_display == "stars"){
      plotting_cells <- cellnote_matrix
    } else if(cellnote_display == "percents"){
      plotting_cells <- apply(100*round(your_data, digits = 3), 2, paste0, "%")
    } else {
      stop("plotting_cells must be one of \"stars\" or \"percents\"")
    }
    if(print_table){
      return(your_data)
    }
    plotting_data <- custom_log(your_data, exp(1), F)
    gplots::heatmap.2(as.matrix(plotting_data),
                      trace = "none", Rowv = F,
                      dendrogram = "none", Colv = FALSE,
                      margins = c(15,6),rowsep = 1:nrow(your_data),
                      colsep = 1:ncol(your_data),
                      keysize = 0.9, sepcolor ="black",
                      cellnote = plotting_cells,
                      notecol = "black",
                      labRow = "",
                      notecex = cellnote_size,
                      cexCol = text_size/10,
                      col = (rainbow(256, s = 1, v = 1, start = 0, end = 0.75, alpha = 1)))

  }else{
    stop("plot_mode must be one of \"line\", \"bar\", or \"heatmap\"")
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
