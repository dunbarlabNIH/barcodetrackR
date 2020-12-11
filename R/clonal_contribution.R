#'@importFrom rlang %||%
#'@importFrom magrittr %>%
#'
#'@title Clonal contribution plot
#'
#'@description Bar or line plot of percentage contribution of the top clones from a selected sample or all clones across samples matching the specified filter within the SummarizedExperiment object. Usually used for tracking a cell lineage's top clones over time.
#'
#'@param your_SE A Summarized Experiment object.
#'@param SAMPLENAME_choice The identifying SAMPLENAME from which to obtain the top "n_clones" clones. If NULL, must specify clone_sequences.
#'@param filter_by Name of metadata column to filter by e.g. Lineage
#'@param filter_selection The value of the filter column to display e.g. "T" (within Lineage)
#'@param plot_over The column of metadata that you want to be the x-axis of the plot. e.g. Month
#'@param plot_over_display_choices Choice(s) from the column designated in plot_over that will be used for plotting. Defaults to all.
#'@param clones_sequences The identifying rownames within your_SE for which to plot. If NULL, must specify SAMPLENAME_choice.
#'@param n_clones Numeric. Number of top clones from SAMPLENAME_choice that should be displayed.
#'@param graph_type Choice of "bar" or "line" for how to display the clonal contribution data
#'@param keep_numeric If plot_over is numeric, whether to space the x-axis appropriately according to the numerical values.
#'@param plot_non_selected Plot clones not in clones_sequences or that aren't top clones in SAMPLENAME_choice
#'@param linesize Numeric. Thickness of the lines.
#'@param text_size Numeric. Size of text in plot.
#'@param y_limit Numeric. What the max value of the y scale should be for the "percentages" assay.
#'@param return_table Logical. If set to TRUE, the function will return a dataframe with each sequence that is selected and its percentage contribution to each selected sample rather than a plot.
#'
#'@importFrom stats setNames
#'
#'@return Displays a stacked area line or bar plot (made by ggplot2) of the samples' top clones. Or, if return_table is set to TRUE, returns a dataframe of the percentage abundances in each sample.
#'
#'@examples
#'clonal_contribution(your_data = your_barcoding_SE, graph_type = "bar",  SAMPLENAME_choice = "Granulocytes_3m", filter_by = "celltype", filter_selection = "Granulocytes", plot_over = "Timepoint_weeks", n_clones = 20)
#'@export

clonal_contribution <- function(your_SE,
                                SAMPLENAME_choice,
                                filter_by,
                                filter_selection,
                                plot_over,
                                plot_over_display_choices = NULL,
                                clone_sequences = NULL,
                                n_clones = 10,
                                graph_type = "bar",
                                keep_numeric = TRUE,
                                plot_non_selected = TRUE,
                                linesize = 0.2,
                                text_size = 15,
                                your_title = "",
                                y_limit = NULL,
                                return_table = FALSE){


  # Some basic error checking before running the function
  coldata_names <- colnames(SummarizedExperiment::colData(your_SE))
  if(any(! c(plot_over) %in% coldata_names)){
    stop("plot_over must match a column name in colData(your_SE)")
  }
  if(sum(is.null(SAMPLENAME_choice), is.null(clone_sequences)) != 1){
    stop("please specify only one of SAMPLENAME_choice or clone_sequences")
  }
  if(any(! c(filter_by) %in% coldata_names)){
    stop("filter_by must match a column name in colData(your_SE)")
  }
  if(! filter_selection %in% unique(SummarizedExperiment::colData(your_SE)[[filter_by]])){
    stop("filter_selection must be an element in the colData column specified with filter_by")
  }

  #get appropriate rows to use in plotting
  if(!is.null(SAMPLENAME_choice)){
    selected_sequences <- get_top_clones(your_SE = your_SE, SAMPLENAME_choice = SAMPLENAME_choice, n_clones = n_clones)
  } else if(!is.null(clone_sequences)){
    selected_sequences <- clone_sequences
  } else {
    stop("one of SAMPLENAME_choice or clone_sequences must be not-NULL")
  }


  #select those samples which to plot_over and to filter_by
  temp_subset <- your_SE[,your_SE[[filter_by]] == filter_selection]
  if(is.numeric(SummarizedExperiment::colData(temp_subset)[[plot_over]])){
    plot_over_display_choices <- plot_over_display_choices  %||% sort(unique(SummarizedExperiment::colData(temp_subset)[[plot_over]]))
    plot_over_display_choices <- as.numeric(as.character(plot_over_display_choices))
  } else if (is.factor(SummarizedExperiment::colData(temp_subset)[[plot_over]])) {
    plot_over_display_choices <- plot_over_display_choices %||% levels(SummarizedExperiment::colData(temp_subset)[[plot_over]])
  } else {
    plot_over_display_choices <- plot_over_display_choices %||% unique(SummarizedExperiment::colData(temp_subset)[[plot_over]])
  }
  temp_subset <- temp_subset[,temp_subset[[plot_over]] %in% plot_over_display_choices]
  temp_subset_coldata <- SummarizedExperiment::colData(temp_subset) %>% as.data.frame() %>% dplyr::mutate_if(is.factor,as.character)


  #ensure that filter_by/filter_selection results in a subset of samples that is identified by a unique element in plot_over
  if(length(temp_subset_coldata[[plot_over]]) != length(unique(temp_subset_coldata[[plot_over]]))){
    stop("after subsetting using filter_by/filter_selection, the remaining elements in the plot_over column must be unique")
  }

  #fetch percentages and turn sample_name into the plot_over equivalents
  your_data <- SummarizedExperiment::assays(temp_subset)[["percentages"]]
  your_data <- your_data[rowSums(your_data) > 0,,drop = FALSE]
  plotting_data <- tibble::rownames_to_column(your_data, var = "sequence") %>%
    tidyr::pivot_longer(cols = -sequence, names_to = "sample_name", values_to = "value") %>%
    dplyr::mutate(fill_label = ifelse(sequence %in% selected_sequences, sequence, "other")) %>%
    dplyr::mutate(fill_label = factor(.data$fill_label, levels = c(selected_sequences, "other"))) %>%
    dplyr::left_join(temp_subset_coldata %>% dplyr::rename(sample_name = .data$SAMPLENAME), by = "sample_name") %>%
    dplyr::mutate(sample_name = !!as.name(plot_over))
  if(is.numeric(temp_subset_coldata[[plot_over]]) & keep_numeric){
    plotting_data <- dplyr::mutate(plotting_data, sample_name = factor(.data$sample_name, levels = plot_over_display_choices))
    plotting_data <- dplyr::mutate(plotting_data, sample_name = as.numeric(as.character(.data$sample_name)))
  } else {
    plotting_data <- dplyr::mutate(plotting_data, sample_name = factor(.data$sample_name, levels = plot_over_display_choices))
  }

  #set up appropriate levels in plotting the selected elements and specify the colors
  if(plot_non_selected){
    non_selected_sequences <- plotting_data %>%
      dplyr::arrange(desc(.data$value)) %>%
      dplyr::pull(sequence) %>%
      unique() %>%
      .[!(. %in% selected_sequences)]
    plotting_data$sequence <- factor(plotting_data$sequence, levels = rev(c(selected_sequences, non_selected_sequences)))
    color_vector <- setNames(c(scales::hue_pal()(length(selected_sequences)), "grey"), c(selected_sequences, "other"))
  } else {
    plotting_data <- dplyr::filter(plotting_data, sequence %in% selected_sequences)
    plotting_data$sequence <- factor(plotting_data$sequence, levels = rev(c(selected_sequences)))
    color_vector <- setNames(c(scales::hue_pal()(length(selected_sequences))), selected_sequences)
  }

  if (return_table){
    return(plotting_data[,1:3])
  }

  if (graph_type == "bar"){
    g <- ggplot2::ggplot(plotting_data, ggplot2::aes(x=.data$sample_name, y = .data$value, group = .data$sequence, fill = .data$fill_label))+
      ggplot2::geom_col(colour = "black",  size= linesize)+
      ggplot2::scale_y_continuous(name = "percentages", labels = function(x){paste0(x * 100, "%")}, expand = c(0.01,0))+
      ggplot2::scale_fill_manual("selected_sequences", values = color_vector)+
      ggplot2::theme_classic()+
      ggplot2::ggtitle(your_title)+
      ggplot2::theme(text = ggplot2::element_text(size=text_size),
                     plot.title = ggplot2::element_text(hjust = 0.5), panel.grid = ggplot2::element_blank())

  } else if (graph_type == "line"){
    g <-  ggplot2::ggplot(plotting_data, ggplot2::aes(x=sample_name, y = value, group = sequence, fill = .data$fill_label))+
      ggplot2::geom_area(position = "stack", colour = "black",  size= linesize)+
      ggplot2::scale_y_continuous(name = "percentages", labels = function(x){paste0(x * 100, "%")}, expand = c(0.01,0))+
      ggplot2::scale_fill_manual("selected_sequences", values = color_vector)+
      ggplot2::theme_classic()+
      ggplot2::ggtitle(your_title)+
      ggplot2::theme(text = ggplot2::element_text(size=text_size))
  }

  if(!is.null(y_limit)){
   g <- g + ggplot2::coord_cartesian(ylim = c(0, y_limit))
  }


  if(is.numeric(temp_subset_coldata[[plot_over]]) & keep_numeric){
    g + ggplot2::scale_x_continuous(paste0(plot_over), breaks = plot_over_display_choices, labels = plot_over_display_choices) +
      ggplot2::theme(legend.position="none")
  } else {
    g + ggplot2::scale_x_discrete(paste0(plot_over), breaks = plot_over_display_choices, labels = plot_over_display_choices) +
      ggplot2::theme(legend.position="none")
  }

}
