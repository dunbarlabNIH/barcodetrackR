#' Autocorrelation Plot
#'
#' Gives the pairwise correlation between each sample-sample pair in the data frame in a group-by-group manner. Considers all samples matching the provided "filter_selection" within the "filter_by" column of your metadata. For each unique value of the "plot_over" argument, plots the correlation of that sample with all other samples in that group as a line plot.
#'
#'@param your_SE Summarized Experiment object containing clonal tracking data as created by the barcodetrackR `create_SE` function.
#'@param method_corr Character. One of "pearson", "spearman", or "kendall".
#'@param filter_by Name of metadata column to filter by e.g. Lineage
#'@param filter_selection The value of the filter column to display e.g. "T" (within Lineage)
#'@param plot_over The column of metadata that you want to be the x-axis of the plot. e.g. Month
#'@param plot_over_display_choices Choice(s) from the column designated in plot_over that will be used for plotting. Defaults to all.
#'@param keep_numeric If plot_over is numeric, whether to space the x-axis appropriately according to the numerical values.
#'@param your_title The title for the plot.
#'@param no_negatives Logical. Whether to make negative correlations = 0.
#'@param point_size Numeric. Size of the points.
#'@param line_size Numeric. Thickness of the lines.
#'@param text_size Numeric. Size of text in plot.
#'@param return_table Logical. If set to true, insetad of returning a plot, the function will return a dataframe containing the pairwise correlation values between samples.
#'
#'@return Plots pairwise correlation plot for the samples in your_SE. Or, a dataframe of the pairwise correlation values if return_table is set to TRUE.
#'
#'@importFrom rlang %||%
#'@importFrom magrittr %>%
#'

autocor_plot = function(your_SE,
                        method_corr ="pearson",
                        filter_by,
                        filter_selection,
                        plot_over,
                        plot_over_display_choices = NULL,
                        keep_numeric = FALSE,
                        your_title = "",
                        no_negatives = FALSE,
                        point_size = 1,
                        line_size = 1,
                        text_size = 15,
                        return_table = FALSE) {

  # Some basic error checking before running the function
  coldata_names <- colnames(SummarizedExperiment::colData(your_SE))
  if(any(! c(plot_over) %in% coldata_names)){
    stop("plot_over must match a column name in colData(your_SE)")
  }
  if(any(! c(filter_by) %in% coldata_names)){
    stop("filter_by must match a column name in colData(your_SE)")
  }
  if(! filter_selection %in% unique(SummarizedExperiment::colData(your_SE)[[filter_by]])){
    stop("filter_selection must be an element in the colData column specified with filter_by")
  }

  #select those samples which to plot_over and to filter_by
  temp_subset <- your_SE[,your_SE[[filter_by]] == filter_selection]
  if(is.numeric(SummarizedExperiment::colData(temp_subset)[[plot_over]])){
    plot_over_display_choices <- plot_over_display_choices %||% sort(unique(SummarizedExperiment::colData(temp_subset)[[plot_over]]))
  } else {
    if(is.factor(SummarizedExperiment::colData(temp_subset)[[plot_over]])){
      plot_over_display_choices <- plot_over_display_choices %||% levels(SummarizedExperiment::colData(temp_subset)[[plot_over]])
    } else {
      plot_over_display_choices <- plot_over_display_choices %||% unique(SummarizedExperiment::colData(temp_subset)[[plot_over]])
    }
  }
  temp_subset <- temp_subset[,temp_subset[[plot_over]] %in% plot_over_display_choices]
  temp_subset_coldata <- SummarizedExperiment::colData(temp_subset) %>% as.data.frame() %>% dplyr::mutate_if(is.factor,as.character)
  #ensure that filter_by/filter_selection results in a subset of samples that is identified by a unique element in plot_over
  if(length(temp_subset_coldata[[plot_over]]) != length(unique(temp_subset_coldata[[plot_over]]))){
    stop("after subsetting using filter_by/filter_selection, the remaining elements in the plot_over column must be unique")
  }

  #extracts proportions assay from your_SE
  temp_subset_coldata  %>%
    dplyr::mutate(my_order = !!as.name(plot_over)) %>%
    dplyr::mutate(my_order = factor(.data$my_order, levels = plot_over_display_choices)) %>%
    dplyr::arrange(.data$my_order) %>%
    dplyr::mutate(my_order = as.character(.data$my_order)) -> sorted_temp_subset_coldata
  plotting_data <- SummarizedExperiment::assays(temp_subset)[["proportions"]]
  colnames(plotting_data) <- plyr::mapvalues(colnames(plotting_data),
                                             from = sorted_temp_subset_coldata$SAMPLENAME,
                                             to = sorted_temp_subset_coldata$my_order)
  plotting_data <- plotting_data[,sorted_temp_subset_coldata$my_order]
  plotting_data_longer <- lapply(1:length(sorted_temp_subset_coldata$my_order), function(i){
    lapply(i:length(sorted_temp_subset_coldata$my_order), function(j){
      temp_df <- data.frame(plotting_data[[i]], plotting_data[[j]])
      temp_df <- temp_df[rowSums(temp_df) > 0, ]
      cortest_results <- cor.test(temp_df[[1]], temp_df[[2]], method = method_corr)
      result_df <- data.frame(sample_i = sorted_temp_subset_coldata$my_order[i],
                              sample_j = sorted_temp_subset_coldata$my_order[j],
                              correlation_value = cortest_results$estimate)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .) %>%
    dplyr::rename(grouping_sample = .data$sample_i, target_sample = .data$sample_j)

  if(is.numeric(plot_over_display_choices) & keep_numeric){
    plotting_data_longer <- dplyr::mutate(plotting_data_longer, target_sample = as.numeric(.data$target_sample))
    plotting_data_longer <- dplyr::mutate(plotting_data_longer, grouping_sample = factor(.data$grouping_sample, levels = plot_over_display_choices))
  } else {
    plotting_data_longer <- dplyr::mutate(plotting_data_longer, target_sample = factor(.data$target_sample, levels = plot_over_display_choices))
    plotting_data_longer <- dplyr::mutate(plotting_data_longer, grouping_sample = factor(.data$grouping_sample, levels = plot_over_display_choices))
  }

  if(no_negatives){
    plotting_data_longer %>%
      dplyr::mutate(correlation_value = ifelse(.data$correlation_value < 0, 0, .data$correlation_value)) -> plotting_data_longer
  }

  if (return_table){
    return(plotting_data_longer)
  }

  gg_autocorplot <- ggplot2::ggplot(plotting_data_longer, ggplot2::aes(x = .data$target_sample, y = .data$correlation_value, group = .data$grouping_sample, color = .data$grouping_sample)) +
    ggplot2::geom_line(size=line_size)+
    ggplot2::geom_point(size=point_size)+
    ggplot2::scale_y_continuous(name = "correlation")+
    ggplot2::theme_classic()+
    ggplot2::scale_color_discrete("sample")+
    ggplot2::ggtitle(your_title)+
    ggplot2::theme(text = ggplot2::element_text(size=text_size))

  if(no_negatives){
    gg_autocorplot <- gg_autocorplot + ggplot2::coord_cartesian(ylim = c(0, 1))
  } else {
    my_min <- min(c(plotting_data_longer$correlation_value, 0))
    gg_autocorplot <- gg_autocorplot + ggplot2::coord_cartesian(ylim = c(my_min,1))
  }

  if(is.numeric(plot_over_display_choices) & keep_numeric){
    gg_autocorplot <-
      gg_autocorplot +
      ggplot2::scale_x_continuous(paste0(plot_over), breaks = plot_over_display_choices, labels = plot_over_display_choices)
  } else {
    gg_autocorplot <-
      gg_autocorplot +
      ggplot2::scale_x_discrete(paste0(plot_over), breaks = plot_over_display_choices, labels = plot_over_display_choices)
  }

  return(gg_autocorplot)


}






