#' LibID Contribution
#'
#' Given a chosen sample, will track the top clones from that samples over the other samples included in the dataset.
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param your_data A data frame of one column, with rownames being the barcodes. Used to get the top clones.
#'@param months Numeric vector. Vector corresponding to the months of each sample.
#'@param n_clones Numeric. Number of top clones from each month that should be displayed.
#'@param linesize Numeric. Thickness of the lines.
#'@param text_size Numeric. Size of text in plot.
#'
#'@param y_limit Numeric. Limit of y axis in plot.
#'@param plot_theme Character. One of "BW", "classic", or "original".
#'@param your_title Your title for the plot.
#'@return Displays a stacked area plot (made by ggplot2) of the samples' top clones.
#'@examples
#'topclonestracker(your_data = zh33_Tlineage, top_clones_choice = zh33T1m, months = c(1,2,3,4,5,6.5, 9.5, 12), y_limit = 50)
#'@export

libID_contribution <- function(your_data, months = c(1:ncol(your_data)), separate_libraries = FALSE, n_clones = 5,
                               linesize = 0.5, text_size = 15, y_limit = 100, plot_theme = "classic",
                               your_title = ""){
  libIDs <- unique(substring(rownames(your_data),1,6))
  colnames(your_data) <- months
  your_data <- as.data.frame(100*prop.table(as.matrix(your_data), margin = 2))
  your_data_list <- list()
  my_top_colors <- c("steelblue2","indianred3", "gold1")
  my_other_colors <- c("lightcyan1", "pink1", "light yellow")
  for(i in 1:length(libIDs)){
    your_data_list[[i]] <- your_data[substring(rownames(your_data),1,6) == libIDs[[i]],]
    your_data_top_clones <- your_data_list[[i]][barcodetrackR::gettopindices(your_data_list[[i]], top = n_clones),]
    your_data_other_clones <- your_data_list[[i]][-barcodetrackR::gettopindices(your_data_list[[i]], top = n_clones),]
    your_data_other_clones <- colSums(your_data_other_clones)
    your_data_list[[i]] <- rbind(your_data_top_clones, your_data_other_clones)
    rownames(your_data_list[[i]])[nrow(your_data_list[[i]])] <- paste0("OTHERS ", libIDs[i])
    your_data_list[[i]] <- reshape2::melt(as.matrix(your_data_list[[i]]))
    colnames(your_data_list[[i]]) <- c("BARCODE", "MONTH", "PROP")
    your_data_list[[i]]$LIBRARY <- ifelse(your_data_list[[i]]$BARCODE == paste0("OTHERS ", libIDs[i]), paste0("OTHERS ", libIDs[i]), paste0("TOP ",substring(your_data_list[[i]]$BARCODE,1,6)))
  }

  myColors <- c("indianred3", "pink1", "gold1", "light yellow", "steelblue2", "lightcyan1")
  
  your_data <- do.call(rbind, your_data_list)
  
  your_data$LIBRARY <- as.factor(your_data$LIBRARY)
  
  your_data$LIBRARY <- factor(your_data$LIBRARY,levels(your_data$LIBRARY)[c(5,2,6,3,4,1)])
  
  names(myColors) <- levels(your_data$LIBRARY)
  
  if(separate_libraries) {
    
    your_data_list <- lapply(1:length(libIDs), function(i) {
      your_data[substring(your_data$LIBRARY,nchar(levels(your_data$LIBRARY)[your_data$LIBRARY]) - 5 ,nchar(levels(your_data$LIBRARY)[your_data$LIBRARY])) == libIDs[i],]
    })
  
    plot_list <- lapply(1:length(libIDs), function(i){

      my_Breaks = pretty(0:ceiling(max(aggregate(your_data_list[[i]]$PROP, by = list(your_data_list[[i]]$MONTH), FUN = sum)[,2])))
      ggplot2::ggplot(your_data_list[[i]], ggplot2::aes(x=MONTH, y = PROP, group = BARCODE, fill = LIBRARY))+
        ggplot2::geom_area(position = "stack", linetype = 1, size = linesize, colour = "black")+
        ggplot2::scale_fill_manual(values = myColors)+
        ggplot2::scale_x_continuous(breaks = seq(1, months[length(months)], by = 1), name = "Month", expand = c(0,0))+
        ggplot2::scale_y_continuous(breaks = my_Breaks,name = "Hematopoietic Contribution %", expand  = c(0,0))+
        ggplot2::theme(text = ggplot2::element_text(size=text_size),
                       axis.ticks.length = ggplot2::unit(.75, "cm"),
                       plot.margin = ggplot2::unit(c(1,1,1.5,1.2),"cm"),
                       axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 20)),
                       axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 20)),
                       panel.background = ggplot2::element_rect(fill = "white", colour = "black"),
                       panel.grid.major = ggplot2::element_line(colour = "black"))+
        ggplot2::ggtitle(paste0(your_title, "\n"))+
        ggplot2::coord_cartesian(ylim = c(0, max(my_Breaks)))
    })
    
    
    do.call(gridExtra::grid.arrange, plot_list)
    
    
  } else {
    ggplot2::ggplot(your_data, ggplot2::aes(x=MONTH, y = PROP, group = BARCODE, fill = LIBRARY))+
      ggplot2::geom_area(position = "stack", linetype = 1, size = linesize, colour = "black")+
      ggplot2::scale_fill_manual(values = myColors)+
      ggplot2::scale_x_continuous(breaks = seq(1, months[length(months)], by = 1), name = "Month", expand = c(0,0))+
      ggplot2::scale_y_continuous(name = "Hematopoietic Contribution %", expand = c(0,0))+
      ggplot2::theme(text = ggplot2::element_text(size=text_size),
                     axis.ticks.length = ggplot2::unit(.75, "cm"),
                     plot.margin = ggplot2::unit(c(1,1,1.5,1.2),"cm"),
                     axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 20)),
                     axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 20)))+
      ggplot2::ggtitle(paste0(your_title, "\n"))+
      ggplot2::coord_cartesian(ylim = c(0,y_limit))
    
    
  }
  
  
  
}
