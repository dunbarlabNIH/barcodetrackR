#' Bias line plot
#'
#' Given a summarized experiment and a specified factor to compare bias between "split_bias_on", shows the value of that bias plotted against another specified factor "split_bias_over" where each clone is represented by a line shaded by its overall abundance in the two samples being compared.
# '
#'@param your_SE SummarizedExperiment of barcode data and associated metadata
#'@param split_bias_on The column of metadata corresponding to cell types (or other factor to be compared.)
#'@param bias_1 The first cell type (or other factor) to be compared. Must be a possible value of the split_bias_on column of your metadata. Will be on the UPPER side of the line plot
#'@param bias_2 The second cell type (or other factor) to be compared. Must be a possible value of the split_bias_on column of your metadata. Will be on the LOWER side of the line plot
#'@param split_bias_over The column of metadata to plot by. If numeric, y axis will be in increasing order. If categorical, it will follow order of metadata.
#'@param bias_over Choice(s) from the column designated in `split_bias_over` that will be used for plotting. Defaults to all.
#'@param remove_unique Logical. If set to true, only clones present in both samples will be considered.
#'@param text_size Numeric. The size of the text in the plot.
#'@param keep_numeric Logical. Whether to keep the numeric spacing within split_bias_over or switch to discrete x scale.
#'@param return_table Logical. If set to TRUE, rather than returnign a plot, the function will return a dataframe containing for each barcode sequence and each point of comparison: the bias value, the added proportion between the two factors at that point (cumul_sum), and the maximum cumul_sum (peak_abundance) of that barcode sequence at any point of comparison.
#'
#'@return Bias line plot for two lineages over time. Or if return_table is set to TRUE, a dataframe containing the bias values for each barcode sequence between the two samples at all points of comparison.
#'
#'@importFrom rlang %||%
#'
#'@examples
#'bias_lineplot(your_se = SE, split_bias_on = "selection_type", bias_1 = "B", bias_2 = "T", split_bias_over = "Timepoint")
#'@export
bias_lineplot <- function(your_SE,
                          split_bias_on,
                          bias_1,
                          bias_2,
                          split_bias_over,
                          bias_over = NULL,
                          remove_unique = FALSE,
                          text_size = 16,
                          keep_numeric = TRUE,
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
  } else {
    bias_over <- bias_over %||% levels(SummarizedExperiment::colData(your_SE)[[split_bias_over]])
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


  lapply(1:length(bias_over), function(i){
    loop_coldata <- temp_subset_coldata %>% dplyr::filter(!!as.name(split_bias_over) %in% bias_over[i])
    loop_bias_1 <- dplyr::filter(loop_coldata, !!as.name(split_bias_on) == bias_1) %>% dplyr::pull("SAMPLENAME") %>% as.character()
    loop_bias_2 <- dplyr::filter(loop_coldata, !!as.name(split_bias_on) == bias_2) %>% dplyr::pull("SAMPLENAME") %>% as.character()
    temp_your_data <- your_data[,c(loop_bias_1, loop_bias_2)]
    temp_your_data <- temp_your_data[rowSums(temp_your_data) > 0,]
    if(remove_unique){
      temp_your_data <- temp_your_data[rowSums(temp_your_data > 0) == 2,]
    }
    temp_bias <- log2((temp_your_data[,loop_bias_1] + 1)/(temp_your_data[,loop_bias_2]+1))
    temp_cumsum <- rowSums(prop.table(as.matrix(temp_your_data), margin = 2))
    return(tibble::tibble(barcode = rownames(temp_your_data), plot_over = bias_over[i], bias = temp_bias, cumul_sum = temp_cumsum))
  }) %>%
    do.call(rbind, .) %>%
    dplyr::group_by(barcode) %>%
    dplyr::mutate(peak_abundance = max(cumul_sum)) %>%
    dplyr::arrange(peak_abundance) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(barcode = factor(barcode, levels = unique(barcode)))-> plotting_data

  if(!keep_numeric & is.numeric(bias_over)){
    plotting_data$plot_over <- factor(plotting_data$plot_over, levels = bias_over)
  }

  if (return_table){
    return(plotting_data)
  }  
  
  g <- ggplot2::ggplot(plotting_data, ggplot2::aes(x = plot_over, y = bias, color = cumul_sum, group = barcode)) +
    ggplot2::geom_line()+
    ggplot2::geom_point()+
    ggplot2::scale_color_gradient(name = "Added Proportions", low = "gray85", high = "black")+
    ggplot2::scale_y_continuous(name = paste0("log bias: log2(", bias_1, "/", bias_2, ")")) +
    ggplot2::theme_classic()+
    ggplot2::theme(text = ggplot2::element_text(size = text_size),
                   axis.text.x = ggplot2::element_text(vjust = 0.5, hjust = 1, angle = 90),
                   legend.position = "right")

  if(keep_numeric & is.numeric(bias_over)){
    g <- g + ggplot2::scale_x_continuous(name = split_bias_over, breaks = bias_over)
  } else {
    g <- g + ggplot2::scale_x_discrete(name = split_bias_over)
  }
  g
}
