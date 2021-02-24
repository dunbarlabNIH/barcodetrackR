#'@importFrom rlang %||%
#'@importFrom magrittr %>%
#'
#'@title Clonal contribution plot
#'
#'@description Bar or line plot of percentage contribution of the top clones from a selected sample or all clones across samples matching the specified filter within the SummarizedExperiment object. Usually used for tracking a cell lineage's top clones over time.
#'
#'@param your_SE A Summarized Experiment object.
#'@param SAMPLENAME_choice The identifying SAMPLENAME from which to obtain the top "n_clones" clones to color. If NULL and clone_sequences is NULL, all clones will be shown as gray.
#'@param filter_by Name of metadata column to filter by e.g. Lineage
#'@param filter_selection The value of the filter column to display e.g. "T" (within Lineage)
#'@param plot_over The column of metadata that you want to be the x-axis of the plot. e.g. Month. For numeric metadata, the x-axis will be ordered in ascending fashion. For categorical metadata, the sample order will be followed.
#'@param plot_over_display_choices Choice(s) from the column designated in plot_over that will be used for plotting. Defaults to all.
#'@param clone_sequences The identifying rownames within your_SE for which to plot. SAMPLENAME_choice should be set to NULL or not specified if clone_sequences is specified.
#'@param n_clones Numeric. Number of top clones from SAMPLENAME_choice that should be assigned a unique color.
#'@param graph_type Choice of "bar" or "line" for how to display the clonal contribution data
#'@param keep_numeric If plot_over is numeric, whether to space the x-axis appropriately according to the numerical values.
#'@param plot_non_selected Plot clones NOT found within the top clones in SAMPLENAME_choice or the specified clones passed to clone_sequences. These clones are colored gray. If both SAMPLENAME_choice and clone_sequences are NULL, this argument must be set to TRUE. Otherwise, there will be no data to show.
#'@param linesize Numeric. Thickness of the lines.
#'@param text_size Numeric. Size of text in plot.
#'@param your_title Title string for your plot.
#'@param y_limit Numeric. What the max value of the y scale should be for the "proportions" assay.
#'@param return_table Logical. If set to TRUE, the function will return a dataframe with each sequence that is selected and its percentage contribution to each selected sample rather than a plot.
#'
#'@importFrom stats setNames
#'
#'@return Displays a stacked area line or bar plot (made by ggplot2) of the samples' top clones. Or, if return_table is set to TRUE, returns a dataframe of the percentage abundances in each sample.
#'
#'@examples
#'clonal_contribution(your_SE = wu_subset, graph_type = "bar",
#'                    SAMPLENAME_choice = "ZJ31_20m_T",
#'                    filter_by = "celltype", filter_selection = "T",
#'                    plot_over = "months", n_clones = 10)
#'@export

clonal_contribution <- function(your_SE,
                                SAMPLENAME_choice = NULL,
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


  # Error checking
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

  if(sum(is.null(SAMPLENAME_choice), is.null(clone_sequences)) == 0){
    stop("please specify only ONE of SAMPLENAME_choice or clone_sequences.")
  }

  if(sum(is.null(SAMPLENAME_choice), is.null(clone_sequences)) == 2 & plot_non_selected == FALSE){
    stop("If neither SAMPLENAME_choice nor clone_sequences are specified, plot_non_selected must be TRUE. Otherwise, there is no data to plot.")
  }

  #get appropriate rows to use in plotting
  if(!is.null(SAMPLENAME_choice)){
    selected_sequences <- get_top_clones(your_SE = your_SE, SAMPLENAME_choice = SAMPLENAME_choice, n_clones = n_clones)
  } else if(!is.null(clone_sequences)){
    selected_sequences <- clone_sequences
  } else {
    selected_sequences <- c("Filler text")
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


  # #ensure that filter_by/filter_selection results in a subset of samples that is identified by a unique element in plot_over
  # if(length(temp_subset_coldata[[plot_over]]) != length(unique(temp_subset_coldata[[plot_over]]))){
  #   stop("after subsetting using filter_by/filter_selection, the remaining elements in the plot_over column must be unique")
  # }

  #fetch proportions
  your_data <- SummarizedExperiment::assays(temp_subset)[["proportions"]]
  your_data <- your_data[rowSums(your_data) > 0,,drop = FALSE]

  # If there are duplicate samples within the chosen plot_over argument for the x-axis, then combine but warn the user
  if(length(temp_subset_coldata[[plot_over]]) != length(unique(temp_subset_coldata[[plot_over]]))){
    duplicated_plot_over <- temp_subset_coldata[[plot_over]][duplicated(temp_subset_coldata[[plot_over]])]
    cat("Duplicate samples with the same value of the plot_over variable:",plot_over, "\n")
    for (i in seq_along(duplicated_plot_over)){
      duplicated_samplenames <- temp_subset_coldata[["SAMPLENAME"]][temp_subset_coldata[[plot_over]] == duplicated_plot_over[i]]
      cat(plot_over, "value =",duplicated_plot_over[i], "; Duplicate sample names =", duplicated_samplenames, "\n")
      your_data[,duplicated_samplenames[1]] <- rowMeans(your_data[,duplicated_samplenames])
      your_data[,duplicated_samplenames[-1]] <- NULL
      temp_subset_coldata <- temp_subset_coldata[temp_subset_coldata[["SAMPLENAME"]] %in% duplicated_samplenames[-1] == FALSE,]
    }
    cat("Barcode proportions have been averaged for duplicate samples.\nTo see the duplicate samples as separate replicates, create a categorical column of metadata with your desired plot_over values but _repX appended to the replicate samples. \n")
  }

  # turn sample_name into the plot_over equivalents
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
    return(plotting_data[,seq_len(3)])
  }

  if (graph_type == "bar"){
    g <- ggplot2::ggplot(plotting_data, ggplot2::aes(x=.data$sample_name, y = .data$value, group = .data$sequence, fill = .data$fill_label))+
      ggplot2::geom_col(colour = "black",  size= linesize)+
      ggplot2::scale_y_continuous(name = "proportions", labels = function(x){paste0(x * 100, "%")}, expand = c(0.01,0))+
      ggplot2::scale_fill_manual("selected_sequences", values = color_vector)+
      ggplot2::theme_classic()+
      ggplot2::ggtitle(your_title)+
      ggplot2::theme(text = ggplot2::element_text(size=text_size),
                     plot.title = ggplot2::element_text(hjust = 0.5), panel.grid = ggplot2::element_blank())

  } else if (graph_type == "line"){
    g <-  ggplot2::ggplot(plotting_data, ggplot2::aes(x=.data$sample_name, y = .data$value, group = sequence, fill = .data$fill_label))+
      ggplot2::geom_area(position = "stack", colour = "black",  size= linesize)+
      ggplot2::scale_y_continuous(name = "proportions", labels = function(x){paste0(x * 100, "%")}, expand = c(0.01,0))+
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
