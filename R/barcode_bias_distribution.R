barcode_bias_distribution <- function(your_data, breaks, CT, TP) {
  your_data <- as.data.frame(prop.table(as.matrix(your_data), margin = 2))
  your_data[is.na(your_data)] <- 0
  arranger <- unlist(lapply(1:TP, function(i){seq(from = i, to = TP*CT, by = TP)}))
  your_data <- your_data[,arranger]
  listy <- list()
  for(i in 1:TP){
    current_timepoint <- your_data[,c((CT*(i-1)+1):(CT*(i-1)+5))]
    current_timepoint <- current_timepoint[current_timepoint[,ncol(current_timepoint)] != 0, ]
    x <- lapply(breaks, function(i){
      current_timepoint_biased_barcodes <- rownames(current_timepoint[apply(current_timepoint[,1:ncol(current_timepoint)-1], 1, max)*i < current_timepoint[,ncol(current_timepoint)],])
      current_timepoint_subset <- current_timepoint[current_timepoint_biased_barcodes,CT,drop = FALSE]
      current_timepoint <<- current_timepoint[!(rownames(current_timepoint) %in% current_timepoint_biased_barcodes),]
      current_timepoint_subset$BARCODES <- rownames(current_timepoint_subset)
      current_timepoint_subset$SAMPLE <- colnames(current_timepoint_subset)[1]
      current_timepoint_subset$FOLD <- paste0(i, "X")
      rownames(current_timepoint_subset) <- NULL
      colnames(current_timepoint_subset)[1] <- "PROP"
      return(current_timepoint_subset)
    })
    listy[[i]] <- do.call(rbind, x)
    listy[[i]]$TIMEPOINT <- paste0("TP#",i)
    print(head(listy[[i]]))
  }
  plotting_data <- do.call(rbind, listy)
  plotting_data$FOLD <- factor(plotting_data$FOLD, levels = paste0(breaks, "X"))
  plotting_data <- plotting_data[order(plotting_data$PROP),]
  ggplot2::ggplot(plotting_data, ggplot2::aes(x = FOLD, y = PROP))+
    ggplot2::geom_bar(ggplot2::aes(y = PROP), stat = "identity",  colour = "black", fill = "white")+
    ggplot2::facet_grid(~ TIMEPOINT)+
    ggplot2::scale_x_discrete(labels = paste0(">", breaks, "X"))+
    ggplot2::scale_y_continuous(expand = c(0,0), breaks = seq(from = 0, to = 1, by = 0.2))+
    ggplot2::coord_cartesian(ylim = c(0,1))+
    ggplot2::guides(color = FALSE, fill = FALSE)+
    ggplot2::theme(text = ggplot2::element_text(size=20),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(colour = "black"),
                   axis.ticks.length = ggplot2::unit(1, "cm"),
                   axis.ticks = ggplot2::element_line(size = 1),
                   panel.background = ggplot2::element_rect(fill = "white", colour = "black"),
                   #panel.grid.major = ggplot2::element_line(colour = "lightgrey"),
                   axis.text.x = ggplot2::element_text(angle = 90))
  
  
  
  
  
  
}