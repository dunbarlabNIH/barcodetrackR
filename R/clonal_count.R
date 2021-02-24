#' Clonal count plot
#'
#' A line plot that tracks the total number of clones or the cumulative number of clones from selected samples of the SummarizedExperiment object plotted over a specified variable.
#'
#'@param your_SE Summarized Experiment object containing clonal tracking data as created by the barcodetrackR `create_SE` function.
#'@param percent_threshold Numeric. The percent threshold for which to count barcodes as present or not present. Set to 0 by default.
#'@param group_by The column of metadata you want to group by e.g. cell_type.
#'@param group_by_choices Choice(s) from the column designated in group_by that will be used for plotting. Defaults to all if left as NULL.
#'@param plot_over The column of metadata that you want to be the x-axis of the plot. e.g. timepoint
#'@param plot_over_display_choices Choice(s) from the column designated in plot_over that will be used for plotting. Defaults to all if left as NULL.
#'@param keep_numeric If plot_over is numeric, whether to space the x-axis appropriately according to the numerical values.
#'@param cumulative Logical. If TRUE, will plot cumulative counts over the `plot_over` argument rather than unique counts per sample (the default, which is FALSE).
#'@param point_size Numeric. Size of points.
#'@param line_size Numeric. Size of lines.
#'@param text_size Numeric. Size of text in plot.
#'@param your_title The title for the plot.
#'@param return_table Logical. If set to true, rather than returning a plot, the function will return the clonal count or cumulative count of each sample in a dataframe.
#'
#'@return Outputs plot of a diversity measure tracked for groups over a factor. Or if return_table is set to TRUE, a dataframe of the number of clones (or cumulative clones) for each sample.
#'
#'@importFrom rlang %||%
#'@importFrom magrittr %>%
#'@importFrom tidyr pivot_longer
#'@import tibble
#'
#'@examples
#'clonal_count(your_SE = wu_subset, cumulative = FALSE, plot_over = "months", group_by = "celltype")
#'
#'@export
clonal_count <- function(your_SE,
                         percent_threshold = 0,
                         plot_over,
                         plot_over_display_choices = NULL,
                         keep_numeric = TRUE,
                         group_by,
                         group_by_choices = NULL,
                         cumulative = FALSE,
                         point_size = 3,
                         line_size = 2,
                         text_size = 12,
                         your_title = NULL,
                         return_table = FALSE) {

  # Some basic error checking before running the function
  coldata_names <- colnames(SummarizedExperiment::colData(your_SE))

  if(!(plot_over %in% coldata_names)){
    stop("plot_over must match a column name in colData(your_SE)")
  }
  if(!(group_by %in% coldata_names)){
    stop("group_by must match a column name in colData(your_SE)")
  }

  if(is.numeric(SummarizedExperiment::colData(your_SE)[[plot_over]])){
    plot_over_display_choices <- plot_over_display_choices  %||% sort(unique(SummarizedExperiment::colData(your_SE)[[plot_over]]))
    plot_over_display_choices <- as.numeric(as.character(plot_over_display_choices))
  } else if (is.factor(SummarizedExperiment::colData(your_SE)[[plot_over]])) {
    plot_over_display_choices <- plot_over_display_choices %||% factor(SummarizedExperiment::colData(your_SE)[[plot_over]],, levels = levels(SummarizedExperiment::colData(your_SE)[[plot_over]]))
  } else {
    plot_over_display_choices <- plot_over_display_choices %||% factor(SummarizedExperiment::colData(your_SE)[[plot_over]], levels = unique(SummarizedExperiment::colData(your_SE)[[plot_over]]))
  }

  group_by_choices <- group_by_choices %||% levels(as.factor(SummarizedExperiment::colData(your_SE)[[group_by]]))

  # More error handling
  if(!all(plot_over_display_choices %in% SummarizedExperiment::colData(your_SE)[,plot_over])){
    stop("All elements of plot_over_display_choices must match cvalues in plot_over column")
  }
  if(!all(group_by_choices %in% SummarizedExperiment::colData(your_SE)[,group_by])){
    stop("All elements of group_by_choices must match values in group_by column")
  }

  # extract metadata
  temp_subset <- your_SE[,(your_SE[[plot_over]] %in% plot_over_display_choices)]
  # Keep only the data included in group_by_choices
  temp_subset <- temp_subset[,(temp_subset[[group_by]] %in% group_by_choices)]
  temp_subset_coldata <- SummarizedExperiment::colData(temp_subset) %>% tibble::as_tibble()
  your_data <- SummarizedExperiment::assays(temp_subset)[["proportions"]]
  your_data <- your_data[rowSums(your_data) > 0, ,drop = FALSE]

  #calculate measure for each sample
  if(!cumulative){
    colSums(your_data > percent_threshold) %>%
      tibble::enframe(name = "SAMPLENAME", value = "index") %>%
      dplyr::mutate(index_type = "unique barcode count") -> calculated_index
    ylabel <- "unique barcode count"
  } else {

    # Get count data into tidy format (and remove zeros)
    your_data %>%
      dplyr::mutate(barcode_seq = rownames(your_data)) %>%
      dplyr::select(.data$barcode_seq, everything()) %>%
      tidyr::pivot_longer(cols = seq_len(ncol(your_data))+1,  names_to = "sample", values_to = "count") %>%
      dplyr::filter(count>0) -> tidy_counts

    # Add group by variable into tidy counts
    tidy_counts$group_by <- temp_subset_coldata[match(tidy_counts$sample,temp_subset_coldata$SAMPLENAME),][[group_by]]

    # Add plot over variable into tidy counts
    tidy_counts$plot_over <- temp_subset_coldata[match(tidy_counts$sample,temp_subset_coldata$SAMPLENAME),][[plot_over]]
    # tidy_counts$plot_over <- factor(tidy_counts$plot_over, levels = levels(plot_over_display_choices))

    # Order tidy counts by desired order
    if(is.numeric(SummarizedExperiment::colData(your_SE)[[plot_over]])){
      tidy_counts %>%
        dplyr::arrange(group_by,plot_over) -> tidy_counts_ordered
    } else {
      tidy_counts %>%
        dplyr::arrange(group_by,factor(plot_over, levels = levels(plot_over_display_choices))) -> tidy_counts_ordered
    }

    # Only keep the first occurence of each barcode within each group_by category
    tidy_counts_ordered %>%
      dplyr::group_by(group_by) %>%
      dplyr::distinct(.data$barcode_seq, .keep_all = TRUE) -> tidy_counts_filtered

    # Summarize number of new barcodes and cumulative barcodes at each timepoint

    if(is.numeric(SummarizedExperiment::colData(your_SE)[[plot_over]])){
      tidy_counts_filtered %>%
        dplyr::group_by(.data$group_by,plot_over,.data$sample) %>%
        dplyr::summarise(new_count = dplyr::n(), .groups = 'drop') %>%
        dplyr::group_by(group_by) %>%
        dplyr::mutate(cumulative_count = cumsum(.data$new_count)) -> summarized_data
    } else {
      tidy_counts_filtered %>%
        dplyr::group_by(group_by,factor(plot_over, levels = levels(plot_over_display_choices)),sample) %>%
        dplyr::summarise(new_count = dplyr::n(), .groups = 'drop') %>%
        dplyr::group_by(group_by) %>%
        dplyr::mutate(cumulative_count = cumsum(.data$new_count)) -> summarized_data
      colnames(summarized_data)[2] <- "plot_over"
    }

    # Put into proper structure
    as.data.frame(summarized_data) %>%
      dplyr::select(.data$sample,.data$cumulative_count) %>%
      dplyr::rename(SAMPLENAME = .data$sample, index = .data$cumulative_count) %>%
      dplyr::mutate(index_type = "unique barcodes") %>%
      dplyr::mutate(SAMPLENAME = as.character(.data$SAMPLENAME)) -> calculated_index

    # Fix the fact that certain samples will be dropped if they have 0 new clones
    # Hacky way to do it. Need to go back and make it better later.
    if (length(unique(tidy_counts_filtered$sample)) < length(unique(tidy_counts_ordered$sample))){
      samples_not_found <- unique(tidy_counts_ordered$sample)[unique(tidy_counts_ordered$sample) %in% unique(tidy_counts_filtered$sample) == FALSE]
      calculated_index <- rbind(calculated_index, data.frame(SAMPLENAME = samples_not_found,
                                                             index = rep(0,length(samples_not_found)),
                                                             index_type = rep("unique barcodes", length(samples_not_found))))
      # Reorder
      calculated_index %>%
        dplyr::arrange(factor(.data$SAMPLENAME), levels = unique(tidy_counts_ordered$sample))  -> calculated_index

      # Give the missing samples the same cumulative count as above samples
      for (i in seq_along(nrow(calculated_index))){
        if (calculated_index[i,"SAMPLENAME"] %in% samples_not_found){
          calculated_index[i,"index"] = calculated_index[i-1,"index"]
        }
      }

    }
    ylabel <- "cumulative barcode count"

  }

  # merge measures with colData
  temp_subset_coldata %>%
    dplyr::mutate(SAMPLENAME = as.character(.data$SAMPLENAME)) %>%
    dplyr::left_join(calculated_index, by = "SAMPLENAME") -> plotting_data

  # Make sure plot over is a factor if not numeric or specified to not keep numeric.
  if(is.numeric(temp_subset_coldata[[plot_over]]) & keep_numeric){
  } else if (is.numeric(temp_subset_coldata[[plot_over]]) & keep_numeric == FALSE){
    plotting_data[[plot_over]] <- factor(plotting_data[[plot_over]], levels = unique(plot_over_display_choices))
  } else {
    plotting_data[[plot_over]] <- factor(plotting_data[[plot_over]], levels = levels(plot_over_display_choices))
  }

  if (return_table){
    return(calculated_index)
  }

  plotting_data$x_value <- plotting_data[[plot_over]]
  plotting_data$group_by <- plotting_data[[group_by]]

  # Create ggplot
  g <- ggplot2::ggplot(plotting_data, ggplot2::aes(x = .data$x_value, y = .data$index, group=.data$group_by, colour=.data$group_by)) +
    ggplot2::geom_line(size = line_size)+
    ggplot2::geom_point(size = point_size)+
    ggplot2::labs(x = plot_over, col = group_by, y = ylabel)+
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size=text_size))+
    ggplot2::ggtitle(your_title)

  if(is.numeric(temp_subset_coldata[[plot_over]]) & keep_numeric){
    g + ggplot2::scale_x_continuous(paste0(plot_over), breaks = plot_over_display_choices, labels = plot_over_display_choices)
  } else {
    g + ggplot2::scale_x_discrete(paste0(plot_over), breaks = plot_over_display_choices, labels = plot_over_display_choices)
  }

}


