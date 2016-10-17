barcode_bias_distribution <- function(your_data, breaks = c(10,5,2,1,0), CT, TP, your_title = "", bias_choice = CT) {
  if(length(breaks) != length(unique(breaks))){
    stop("breaks must be unique")
  }
  if(any(breaks[order(-breaks)] != breaks)){
    stop("breaks must be provided from greatest to least")
  }
  your_data <- as.data.frame(prop.table(as.matrix(your_data), margin = 2))
  your_data[is.na(your_data)] <- 0
  listy <- list()
  for(i in 1:TP){
    current_timepoint <- your_data[,c((CT*(i-1)+1):(CT*i)),drop=FALSE]
    current_timepoint <- current_timepoint[current_timepoint[,bias_choice] != 0,]
    x <- lapply(breaks, function(j){
      temp <- current_timepoint[apply(current_timepoint[,-c(bias_choice)], 1, max)*j < current_timepoint[,bias_choice],]
      if(nrow(temp) == 0){
        current_timepoint_biased_barcodes = ""
      } else {
        current_timepoint_biased_barcodes <- rownames(temp)
      }
      current_timepoint_subset <- current_timepoint[current_timepoint_biased_barcodes,bias_choice,drop = FALSE]
      current_timepoint <<- current_timepoint[!(rownames(current_timepoint) %in% current_timepoint_biased_barcodes),]
      current_timepoint_subset$BARCODES <- rownames(current_timepoint_subset)
      current_timepoint_subset$SAMPLE <- colnames(current_timepoint_subset)[1]
      current_timepoint_subset$FOLD <- paste0(j, "X")
      rownames(current_timepoint_subset) <- NULL
      colnames(current_timepoint_subset)[1] <- "FRACTION"
      return(current_timepoint_subset)
    })
    listy[[i]] <- do.call(rbind, x)
    listy[[i]]$TIMEPOINT <- paste0("TP#",i)
  }
  plotting_data <- do.call(rbind, listy)
  plotting_data$FOLD <- factor(plotting_data$FOLD, levels = paste0(breaks, "X"))
  plotting_data <- plotting_data[order(plotting_data$FRACTION),]
  ggplot2::ggplot(plotting_data, ggplot2::aes(x = FOLD, y = FRACTION))+
    ggplot2::geom_bar(ggplot2::aes(y = FRACTION), stat = "identity",  colour = "black", fill = "white")+
    ggplot2::facet_grid(~ TIMEPOINT)+
    ggplot2::scale_x_discrete(labels = paste0(">", breaks, "X"), name = "\nFold Bias\n")+
    ggplot2::scale_y_continuous(expand = c(0,0), breaks = seq(from = 0, to = 1, by = 0.2), labels = function(x){paste0(x*100, "%")})+
    ggplot2::coord_cartesian(ylim = c(0,1))+
    ggplot2::guides(color = FALSE, fill = FALSE)+
    ggplot2::ggtitle(paste0('\n', your_title, '\n'))+
    ggplot2::theme(text = ggplot2::element_text(size=20),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(colour = "black"),
                   axis.ticks.length = ggplot2::unit(1, "cm"),
                   axis.ticks = ggplot2::element_line(size = 1),
                   panel.background = ggplot2::element_rect(fill = "white", colour = "black"),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))






}
