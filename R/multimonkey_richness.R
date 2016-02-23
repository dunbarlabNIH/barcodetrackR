multimonkey_richness_plot <- function(outfile_list,
                                      outfile_names,
                                      outfile_months_list,
                                      celltypes_list,
                                      threshold_list,
                                      combine = FALSE,
                                      richness_type = "unique",
                                      point_size = 5, line_size = 3,
                                      y_lower = 0, y_upper = 9.5){


  for(i in 1:length(outfile_list)){
    outfile_list[[i]] <- barcodetrackR::richness_plot(outfile_list[[i]],
                                                      months = outfile_months_list[[i]],
                                                      celltypes = celltypes_list[[i]],
                                                      thresh = threshold_list[[i]],
                                                      combine = combine,
                                                      richness_type = richness_type,
                                                      show_table = TRUE)
    outfile_list[[i]]$CELLTYPE <- paste(outfile_list[[i]]$CELLTYPE, outfile_names[[i]])

  }
  outfile_list <- do.call(outfile_list, rbind)

  ggplot(outfile_list, aes(x = MONTH, y = BARCODES, group = CELLTYPE, colour = CELLTYPE))+
    geom_line(size = line_size)+
    geom_point(size = point_size)+
    ylab(paste0(richness_type, " Barcodes"))+
    xlab("Month")+
    coord_cartesian(ylim = c(y_lower, y_upper))+
    theme_bw()







}
