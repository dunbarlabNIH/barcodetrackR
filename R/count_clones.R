#' Clonal count plot
#'
#' A line plot that tracks a diversity measure from a selected choice or number of elements in the rows of the SummarizedExperiment object.
#'
#'@param your_SE A Summarized Experiment object.
#'@param group_by The column of metadata you want to group by e.g. cell_type
#'@param plot_over The column of metadata that you want to be the x-axis of the plot. e.g. timepoint
#'@param plot_over_display_choices Choice(s) from the column designated in plot_over that will be used for plotting. Defaults to all if left as NULL.
#'@param index_type Character. One of "count" or "cumulative_count". 
#'@param point_size Numeric. Size of points.
#'@param line_size Numeric. Size of lines.
#'@param text_size Numeric. Size of text in plot.
#'@param your_title The title for the plot.
#'
#'@return Outputs plot of a diversity measure tracked for groups over a factor.
#'
#'@importFrom rlang %||%
#'@importFrom magrittr %>%
#'@importFrom tidyr pivot_longer
#'
#'@examples
#'count_clones(your_data = wu_SE, index_type = "shannon", plot_by = timepoint, group_by = cell_type)
#'
#'@export
count_clones <- function(your_SE,
                         plot_over,
                         plot_over_display_choices = NULL,
                         group_by,
                         index_type = "shannon",
                         point_size = 3,
                         line_size = 2,
                         text_size = 12,
                         your_title = NULL) {

  # Some basic error checking before running the function
  coldata_names <- colnames(SummarizedExperiment::colData(your_SE))
  if(!(plot_over %in% coldata_names)){
    stop("plot_over must match a column name in colData(your_SE)")
  }
  if(!(group_by %in% coldata_names)){
    stop("group_by must match a column name in colData(your_SE)")
  }
  if(is.numeric(SummarizedExperiment::colData(your_SE)[[plot_over]])){
    plot_over_display_choices <- plot_over_display_choices %||% sort(unique(SummarizedExperiment::colData(your_SE)[[plot_over]]))
  } else {
    plot_over_display_choices <- plot_over_display_choices %||% levels(SummarizedExperiment::colData(your_SE)[[plot_over]])
  }

  # extract metadata
  temp_subset <- your_SE[,(your_SE[[plot_over]] %in% plot_over_display_choices)]
  temp_subset_coldata <- SummarizedExperiment::colData(temp_subset) %>% tibble::as_tibble()
  your_data <- SummarizedExperiment::assays(temp_subset)[["percentages"]]
  your_data <- your_data[rowSums(your_data) > 0, ,drop = FALSE]

  #calculate measure for each sample
  if(index_type == "count"){
    colSums(your_data > 0) %>%
      tibble::enframe(name = "SAMPLENAME", value = "index") %>%
      dplyr::mutate(index_type = index_type) -> calculated_index
  } else if (index_type == "cumulative_count"){
    
    # Get count data into tidy format (and remove zeros)
    your_data %>%
      dplyr::mutate(barcode_seq = rownames(your_data)) %>% 
      dplyr::select(barcode_seq, everything()) %>% 
      tidyr::pivot_longer(cols = 1:ncol(your_data)+1,  names_to = "sample", values_to = "count") %>% 
      filter(count>0) -> tidy_counts
    
    # Add group by variable into tidy counts
    tidy_counts$group_by <- temp_subset_coldata[match(tidy_counts$sample,temp_subset_coldata$SAMPLENAME),][[group_by]]
    
    # Add plot over variable into tidy counts
    tidy_counts$plot_over <- temp_subset_coldata[match(tidy_counts$sample,temp_subset_coldata$SAMPLENAME),][[plot_over]]
    
    # Order tidy counts by desired order
    tidy_counts %>%
      arrange(factor(plot_over, levels = plot_over_display_choices),group_by) -> tidy_counts_ordered
    
    # Only keep the first occurence of each barcode within each group_by category
    tidy_counts_ordered %>%
      dplyr::group_by(group_by) %>% 
      dplyr::distinct(barcode_seq, .keep_all = TRUE) -> tidy_counts_filtered
    
    # Summarize number of new barcodes and cumulative barcodes at each timepoint
    tidy_counts_filtered %>% 
      dplyr::group_by(group_by,plot_over,sample) %>% 
      dplyr::summarise(new_count = n(), .groups = 'drop') %>% 
      dplyr::group_by(group_by) %>% 
      dplyr::mutate(cumulative_count = cumsum(new_count)) -> summarized_data
    
    # Put into proper structure
    as.data.frame(summarized_data) %>% 
      dplyr::select(sample,cumulative_count) %>% 
      dplyr::rename(SAMPLENAME = sample, index = cumulative_count) %>% 
      dplyr::mutate(index_type = index_type) -> calculated_index
    
  } else {
    stop("index_type must be one of \"count\", \"cumulative_count\"")
  }

  # merge measures with colData
  temp_subset_coldata %>%
    dplyr::mutate(SAMPLENAME = as.character(SAMPLENAME)) %>%
    dplyr::left_join(calculated_index, by = "SAMPLENAME") -> plotting_data

  plotting_data[[plot_over]] <- factor(plotting_data[[plot_over]], levels = plot_over_display_choices)
  plotting_data$x_value <- plotting_data[[plot_over]]
  plotting_data$group_by <- plotting_data[[group_by]]

  # Create ggplot
  ggplot2::ggplot(plotting_data, ggplot2::aes(x = x_value, y = index, group=group_by, colour=group_by)) +
    ggplot2::geom_line(size = line_size)+
    ggplot2::geom_point(size = point_size)+
    ggplot2::labs(x = plot_over, col = group_by, y = gsub("_"," ",index_type))+
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size=text_size))+
    ggplot2::ggtitle(your_title)
}


