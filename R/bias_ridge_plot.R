#' Bias Ridge plot
#'
#' Given a summarized experiment and a specified factor to compare bias between, gives ridge plots which show the density of clones at each value of log bias where log bias is calculated as log((normalized abundance in sample 1 + 1)/(normalized abundance in sample 2 + 1)). If the weighted option is set to TRUE, the density estimator will weight the estimation by the added proportion of the clone between the two samples.
# '
#'@param your_SE Your SummarizedExperiment of barcode data and associated metadata
#'@param split_bias_on The column of metadata corresponding to cell types (or whatever factors you want to compare the bias between).
#'@param bias_1 The first cell type (or other factor) to be compared. Must be a possible value of the split_bias_on column of your metadata. Will be on the RIGHT side of the ridge plot
#'@param bias_2 The second cell type (or other factor) to be compared. Must be a possible value of the split_bias_on column of your metadata. Will be on the LEFT side of the ridge plot
#'@param split_bias_over The column of metadata to plot by. If numeric, y axis will be in increasing order. If categorical, it will follow order of metadata.
#'@param bias_over Choice(s) from the column designated in `split_bias_over` that will be used for plotting. Defaults to all.
#'@param remove_unique If set to true, only clones present in both samples will be considered.
#'@param weighted If true, the density estimation will be weighted by the overall contribution of each barcode to the two samples being compared.
#'@param text_size Numeric. The size of the text in the plot.
#'@param add_dots Logical. Whether or not to add dots underneath the density plots. Dot size is proportion to the added proportion of the clone in the two samples.
#'@param return_table Logical. If true, rather than returning a plot, the function will return a dataframe containing the calculated bias and cumul_sum which contains the added proportion between the two samples, for each barcode sequence across each sample considered.
#'
#'@return Bias plot for two lineages over time. Or a dataframe containing the bias value and added proportion of each barcode if return_table is set to TRUE.
#'
#'@importFrom rlang %||%
#'@importFrom plyr .
#'@import ggridges
#'@import utils
#'
#'@examples
#'bias_ridge_plot(your_SE = wu_subset, split_bias_on = "celltype",
#'                bias_1 = "B", bias_2 = "T", split_bias_over = "months",
#'                add_dots = TRUE)
#'
#'@export
bias_ridge_plot <- function(your_SE,
                       split_bias_on,
                       bias_1,
                       bias_2,
                       split_bias_over,
                       bias_over = NULL,
                       remove_unique = FALSE,
                       weighted = FALSE,
                       text_size = 16,
                       add_dots = FALSE,
                       return_table = FALSE){

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
  } else if (is.factor(SummarizedExperiment::colData(your_SE)[[split_bias_over]])){
    bias_over <- bias_over %||% levels(SummarizedExperiment::colData(your_SE)[[split_bias_over]])
  } else {
    bias_over <- bias_over %||% unique(SummarizedExperiment::colData(your_SE)[[split_bias_over]])
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
  your_data <- SummarizedExperiment::assays(temp_subset)[["normalized"]]
  your_data <- your_data[rowSums(your_data) > 0, ,drop = FALSE]


  lapply(seq_along(bias_over), function(i){
    loop_coldata <- temp_subset_coldata %>% dplyr::filter(!!as.name(split_bias_over) %in% bias_over[i])
    loop_bias_1 <- dplyr::filter(loop_coldata, !!as.name(split_bias_on) == bias_1) %>% dplyr::pull("SAMPLENAME") %>% as.character()
    loop_bias_2 <- dplyr::filter(loop_coldata, !!as.name(split_bias_on) == bias_2) %>% dplyr::pull("SAMPLENAME") %>% as.character()
    temp_your_data <- your_data[,c(loop_bias_1, loop_bias_2)]
    temp_your_data <- temp_your_data[rowSums(temp_your_data) > 0,]
    if(remove_unique){
      temp_your_data <- temp_your_data[rowSums(temp_your_data > 0) == 2,]
    }
    temp_bias <- log2((temp_your_data[,loop_bias_1] + 1)/(temp_your_data[,loop_bias_2]+1))
    temp_cumsum <- rowSums(temp_your_data)
    return(tibble::tibble(barcode = rownames(temp_your_data), plot_over = bias_over[i], bias = temp_bias, cumul_sum = temp_cumsum))
  }) %>%
    do.call(rbind, .) %>%
    dplyr::mutate(plot_over = factor(.data$plot_over, levels = bias_over)) -> plotting_data

  if (return_table){
    return(plotting_data)
  }

  # Weighted ridge plot
  if (weighted){
    g <- ggplot2::ggplot(plotting_data, ggplot2::aes(x = .data$bias, y = .data$plot_over, height = .data$..density.., weight = .data$cumul_sum, fill = .data$plot_over))
  } else {
    g <- ggplot2::ggplot(plotting_data, ggplot2::aes(x = .data$bias, y = .data$plot_over, height = .data$..density.., fill = .data$plot_over))
  }
  if(add_dots){
    g <- g +
      ggplot2::geom_point(data = plotting_data, ggplot2::aes(x = .data$bias, y = .data$plot_over, size = .data$cumul_sum), inherit.aes = FALSE)
  }

  g +
    ggridges::geom_density_ridges(stat = 'density', alpha = 0.5) +
    ggplot2::scale_x_continuous(name = paste0("log bias: log2(", bias_1, "/", bias_2, ")")) +
    ggplot2::theme(text = ggplot2::element_text(size = text_size),
                   panel.grid.major.x = ggplot2::element_line(colour = "grey"),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.background = ggplot2::element_rect(fill = "white"),
                   axis.text.x = ggplot2::element_text(vjust = 0.5, hjust = 1, angle = 90),
                   legend.position = "none") + ggplot2::scale_y_discrete(name = split_bias_over)

}
