#'@title Bias histogram
#'
#'@description Given a summarized experiment, gives histogram of log biases for 2 cell types. Each stacked bar in the histogram represents a clone binned by log bias defined as the log2 of the percentage abundance in the sample specified in "bias_1" divided by the percentage abundance in "bias_2."
#'
#'@param your_SE Your SummarizedExperiment of barcode data and associated metadata
#'@param split_bias_on The column in `colData(your_SE)` from which `bias_1` and `bias_2` will be chosen
#'@param bias_1 The first cell type (or other factor) to be compared. Must be a possible value of the split_bias_on column of your metadata. Will be on the RIGHT side of the histogram.
#'@param bias_2 The second cell type (or other factor) to be compared. Must be a possible value of the split_bias_on column of your metadata. Will be on the LEFT side of the ridge plot
#'@param split_bias_over The column in `colData(your_SE)` that you wish to observe the bias split for. The output will contain a faceted plot: one facet for each value of `split_bias_over` comparing the samples matching `bias_1` and `bias_2` from the `split_bias_on` argument.
#'@param bias_over Choice(s) from the column designated in `split_bias_over` that will be used for plotting. Defaults to all.
#'@param remove_unique If set to true, only clones present in both samples will be considered.
#'@param breaks Numeric. The breaks specified for bins on the x-axis (how biased the clones are towards one factor or the other).
#'@param text_size The size of the text in the plot.
#'@param linesize The linewidth of the stacked bars which represent individual barcodes
#'@param ncols Numeric. Number of columns to plot on using plot_grid from cowplot.
#'@param scale_all_y Logical. Whether or not to plot all plots on the same y axis limits.
#'@param return_table Logical. If set to TRUE, instead of a plot, tbe function will return a list containing a dataframe for each sample-sample log bias combination containing each barcode sequence and its bias between the samples.
#'
#'@return Histogram of log bias for two factors faceted over another set of factors. Or, if return_table is set to TRUE, a list of dataframes containing the log bias data for each bias comparison passed to the function.
#'
#'@importFrom rlang %||%
#'@importFrom magrittr %>%
#'
#'@examples
#'bias_histogram(your_SE = wu_subset, split_bias_on = "celltype",
#'               bias_1 = "B", bias_2 = "T",
#'               split_bias_over = "months", ncols = 2)
#'@export
bias_histogram <- function(your_SE,
                           split_bias_on,
                           bias_1,
                           bias_2,
                           split_bias_over,
                           bias_over = NULL,
                           remove_unique = FALSE,
                           breaks = c(10,2,1,0.5),
                           text_size = 10,
                           linesize = .4,
                           ncols = 1,
                           scale_all_y = TRUE,
                           return_table = FALSE) {

  breaks_labels <- breaks
  breaks <- c(-Inf, sort(unique(c(-breaks, breaks))), Inf)

  # Some basic error checking before running the function
  if(length(breaks) != length(unique(breaks))){
    stop("breaks must be unique")
  }
  # Basic error handling
  coldata_names <- colnames(SummarizedExperiment::colData(your_SE))
  if(any(! c(split_bias_over) %in% coldata_names)){
    stop("split_bias_over must match a column name in colData(your_SE)")
  }
  if(any(! c(split_bias_on) %in% coldata_names)){
    stop("split_bias_on must match a column name in colData(your_SE)")
  }
  if(! bias_1 %in% unique(SummarizedExperiment::colData(your_SE)[[split_bias_on]])){
    stop("bias_1 is not in split_bias_on")
  }
  if(! bias_2 %in% unique(SummarizedExperiment::colData(your_SE)[[split_bias_on]])){
    stop("bias_2 is not in split_bias_on")
  }
  if(! all(bias_over %in% unique(SummarizedExperiment::colData(your_SE)[[split_bias_over]]))){
    stop("An element in bias_over is not in split_bias_over")
  }

  #perform ordering for the split_bias_over element if numeric and set the variable if it was initially NULL
  if(is.numeric(SummarizedExperiment::colData(your_SE)[[split_bias_over]])){
    bias_over <- bias_over %||% sort(unique(SummarizedExperiment::colData(your_SE)[[split_bias_over]]))
  } else {
    bias_over <- bias_over %||% levels(SummarizedExperiment::colData(your_SE)[[split_bias_over]])
  }

  # Check each value of bias_over
  for (i in seq_along(bias_over)){
    num_samples1 <- nrow(SummarizedExperiment::colData(your_SE)[SummarizedExperiment::colData(your_SE)[[split_bias_over]] == bias_over[i] & SummarizedExperiment::colData(your_SE)[[split_bias_on]] == bias_1,])
    num_samples2 <- nrow(SummarizedExperiment::colData(your_SE)[SummarizedExperiment::colData(your_SE)[[split_bias_over]] == bias_over[i] & SummarizedExperiment::colData(your_SE)[[split_bias_on]] == bias_2,])

    # If some values don't have a comparison, let the user know
    if (num_samples1 == 0 | num_samples2 == 0){
      cat("Note: For bias_over variable", split_bias_over, "the element", bias_over[i], "is missing one or both of the comparators: \n")
      if (num_samples1 == 0){
        cat(split_bias_on,bias_1, "not found. \n \n")
      }
      if (num_samples2 == 0){
        cat(split_bias_on,bias_2, "not found. \n \n")
      }
    }

    # If some values have multiple replicates, let the user know.
    if (num_samples1 > 1 | num_samples2 > 1){
      cat("For bias_over variable", split_bias_over, "the element", bias_over[i], "has more than one replicate for one or more of the comparators. \n")
      if (num_samples1 > 1){
        cat(split_bias_on, bias_1, "has", num_samples1, "replicates. \n \n")
      }
      if (num_samples2 > 1){
        cat(split_bias_on, bias_2, "has", num_samples2, "replicates. \n \n")
      }
    }
  }

  # Repeat the loop to throw the error if there are ambiguous samples after printing all helpful info.
  for (i in seq_along(bias_over)){
    num_samples1 <- nrow(SummarizedExperiment::colData(your_SE)[SummarizedExperiment::colData(your_SE)[[split_bias_over]] == bias_over[i] & SummarizedExperiment::colData(your_SE)[[split_bias_on]] == bias_1,])
    num_samples2 <- nrow(SummarizedExperiment::colData(your_SE)[SummarizedExperiment::colData(your_SE)[[split_bias_over]] == bias_over[i] & SummarizedExperiment::colData(your_SE)[[split_bias_on]] == bias_2,])
    if (num_samples1 > 1 | num_samples2 > 1){
      stop("In order to ensure that the function compares the correct samples, please disambiguate the samples by creating a column of the metadata with _repX appended to the desired `split_bias_over` variable.")
    }
  }

  # ensure that the  chosen bias_over only is able to plot elements in which the chosen bias_1 and bias_2 are present at n = 1 each
  # within each split_bias_over choice
  colData(your_SE) %>%
    tibble::as_tibble() %>%
    dplyr::filter(!!as.name(split_bias_over) %in% bias_over) %>%
    dplyr::filter(!!as.name(split_bias_on) %in% c(bias_1, bias_2)) %>%
    dplyr::group_by(!!as.name(split_bias_over)) %>%
    dplyr::filter((dplyr::n() == 2)) %>%
    dplyr::filter(any(!!as.name(split_bias_on) %in% bias_1) & any(!!as.name(split_bias_on) %in% bias_2))  -> temp_subset_coldata
  bias_over_possibilities <- dplyr::group_keys(temp_subset_coldata) %>% dplyr::pull(!!as.name(split_bias_over))
  bias_over <- bias_over[bias_over %in% bias_over_possibilities]

  # extract bc data
  temp_subset <- your_SE[,temp_subset_coldata$SAMPLENAME]
  your_data <- SummarizedExperiment::assays(temp_subset)[["proportions"]]
  your_data <- your_data[rowSums(your_data) > 0, ,drop = FALSE]

  # pre-allocate
  your_data_list <- list()

  plot_list <- lapply(seq_along(bias_over), function(i){
    loop_coldata <- temp_subset_coldata %>% dplyr::filter(!!as.name(split_bias_over) %in% bias_over[i])
    loop_bias_1 <- dplyr::filter(loop_coldata, !!as.name(split_bias_on) == bias_1) %>% dplyr::pull("SAMPLENAME") %>% as.character()
    loop_bias_2 <- dplyr::filter(loop_coldata, !!as.name(split_bias_on) == bias_2) %>% dplyr::pull("SAMPLENAME") %>% as.character()
    temp_your_data <- your_data[,c(loop_bias_1, loop_bias_2)]
    temp_your_data <- temp_your_data[rowSums(temp_your_data) > 0,]
    if(remove_unique){
      temp_your_data <- temp_your_data[rowSums(temp_your_data > 0) == 2,]
    }
    colnames(temp_your_data) <- c("bias_1", "bias_2")
    temp_your_data %>%
      tibble::rownames_to_column(var = "barcode") %>%
      dplyr::mutate(added_proportions = .data$bias_1 + .data$bias_2, bias = .data$bias_1/.data$bias_2) %>%
      dplyr::mutate(log2_bias = log2(.data$bias)) %>%
      dplyr::mutate(log2_bias_cuts = cut(.data$log2_bias, breaks = breaks, include.lowest = TRUE)) -> temp_your_data

    if (return_table){
      me <- temp_your_data
    } else {
    g <- ggplot2::ggplot(temp_your_data[order(temp_your_data$added_proportions),],
                         ggplot2::aes(x = .data$log2_bias_cuts, y = .data$added_proportions))+
      ggplot2::geom_bar(stat = "identity", fill = "white", size = linesize, color = "black")+
      ggplot2::scale_x_discrete(name = paste0("log bias: log2(", bias_1, "/", bias_2, ")"), drop = FALSE)+
      ggplot2::scale_y_continuous(name = "Added Proportions", labels = function(i)(paste0(i*100, "%")))+
      ggplot2::theme(text = ggplot2::element_text(size = text_size),
                     panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_line(colour = "grey"),
                     panel.background = ggplot2::element_rect(fill = "white", colour = "black"),
                     panel.spacing = ggplot2::unit(2, "lines"),
                     axis.text.x = ggplot2::element_text(vjust = 0.5, hjust = 1, angle = 90))+
      ggplot2::labs(title = paste0(split_bias_over,": ", bias_over[i]))+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }


  })

  if (return_table){
    return(plot_list)
  }

  if(scale_all_y){
    plot_max <- max(unlist(lapply(plot_list, function(x){ggplot2::ggplot_build(x)$layout$panel_params[[1]]$y.range[2]})))
    plot_list <- lapply(plot_list, function(x){x + ggplot2::coord_cartesian(ylim = c(0, plot_max))})
  }

  cowplot::plot_grid(plotlist = plot_list, ncol = ncols)
}


