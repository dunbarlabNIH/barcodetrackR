#' MultiMonkey Richness
#' Plots (using ggplot2) a line plot that tracks the barcode counts of cell types across time points and monkeys. Ideally should only be
#' used with the app, as the data set up can be quite cumbersome otherwise.
#'
#'@param outfile_list A list of outfiles as data frames. See note below for organization of each data frame.
#'@param outfile_names A list of names that correspond to the monkey names of your outfiles.
#'@param outfile_months_list A list of the months for each of the corresponding data frames.
#'@param celltypes_list A list of the celltypes for each of the corresponding data frames.
#'@param threshold_list A list of the thresholds for each of the corresponding data frames.
#'@param combine Logical. Whether to combine the monthly barcodes.
#'@param point_size Numeric. Size of points.
#'@param line_size Numeric. Size of lines.
#'@param richness_type Character. One of "unique", "cumulative", or "new".
#'@param y_lower Numeric. Lower bound of y axis.
#'@param y_upper Numeric. Upper bound of y axis.
#'@return Outputs plot of barcode richness tracked over time.
#'@note Each data frame must be organized in columns by month, and then by cell type. For example, when tracking B and T cells
#'over three months, your columns should be in this order: "1m T", "1m B", "2m T", "2m B", "3m T", "3m B". Then
#'celltypes should be c("T", "B"), and months should be c(1,2,3). Repeat this for every monkey.
#'@examples
#'multimonkey_richness_plot(outfile_list = list(zh33, zg66, zh19),
#'                          outfile_names = list("zh33", "zg66", "zh19"),
#'                          celltypes_list = list(c("T", "B", "Gr"), c("T"), c("T, "B", "Gr)),
#'                          threshold_list = list(1000,1000,0),
#'                          combine = FALSE,
#'                          richness_type = "new")
#'
#'
#'@export
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

  outfile_list <- do.call(rbind, outfile_list)

  ggplot2::ggplot(outfile_list, ggplot2::aes(x = MONTH, y = BARCODES, group = CELLTYPE, colour = CELLTYPE))+
    ggplot2::geom_line(size = line_size)+
    ggplot2::geom_point(size = point_size)+
    ggplot2::ylab(paste0(richness_type, " Barcodes"))+
    ggplot2::xlab("Month")+
    ggplot2::coord_cartesian(ylim = c(y_lower, y_upper))+
    ggplot2::theme_bw()







}
