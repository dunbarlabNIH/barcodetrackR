megatopclones <- function(your_data, num_celltypes, n_clones = 10, threshold = 0){
  rownames(your_data) <- c(1:nrow(your_data))
  your_data[your_data < threshold] <- 0
  your_data <- data.frame(prop.table(as.matrix(your_data),margin = 2))
  split_data <- split(data.frame(t(your_data)), rep(1:num_celltypes, each = (ncol(your_data)/num_celltypes)))
  for (i in 1:length(split_data)){
    split_data[[i]] <- t(split_data[[i]])
    rownames(split_data[[i]]) = c(1:nrow(split_data[[i]]))
  }
  top_clones_per_lineage <- lapply(split_data, gettopindices, top = n_clones)
  top_clones_common <- unlist(top_clones_per_lineage)
  names(top_clones_common) <- NULL
  top_clones_common <- unique(top_clones_common[duplicated(top_clones_common)])
  myColors <- rainbow(length(top_clones_common))
  names(myColors) <- top_clones_common
  myColors <- c(myColors, "N" = "grey")
  split_data_topclones <- list()
  for (i in c(1:num_celltypes)){
    split_data_topclones[[i]] <- lapply(split_data, "[", top_clones_per_lineage[[i]],)
    split_data_topclones[[i]] <- lapply(split_data_topclones[[i]], reshape2::melt)
  }
  split_data_topclones <- unlist(split_data_topclones, recursive = FALSE)
  for (i in c(1:length(split_data_topclones))){
    split_data_topclones[[i]]$CLASS <- ifelse((split_data_topclones[[i]]$Var1 %in% top_clones_common), yes = split_data_topclones[[i]]$Var1,
                                              no = "N")
  }
  plot_list <- list()
  for (i in c(1:length(split_data_topclones))){
    plot_list[[i]] <- ggplot2::ggplot(split_data_topclones[[i]], ggplot2::aes(x = Var2, y = value, group = Var1, fill = CLASS))+
      ggplot2::geom_area(position = "stack", linetype = 1, size = 0.05, colour = "black")+
      ggplot2::theme(legend.position = "none")+
      ggplot2::scale_fill_manual(values = myColors)+
      ggplot2::coord_cartesian(ylim = c(0,0.25))
  }

  do.call(gridExtra::grid.arrange, plot_list)

}
