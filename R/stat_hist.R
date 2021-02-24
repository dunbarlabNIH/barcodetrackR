#' Stat histogram
#'
#' Given a summarized experiment, gives a histogram of the acc assay or choice of metadata.
# '
#'@param your_SE Your SummarizedExperiment of barcode data and associated metadata.
#'@param data_choice Either "assay stats" which allows you to view the distribution of values in the `assay_choice` assay, or "metadata stats" which allows you to view the distribution of metadata values in your SummarizedExperiment object.
#'@param assay_choice When data_choice is set to "assay stats", designates which assay will be used.
#'@param metadata_stat When data_choice is set to "metadata stats", The metadata values that will be used.
#'@param group_meta_by When data_choice is set to "metadata stats", facet the histogram using this column of metadata. If NULL, no grouping or faceting applied
#'@param scale_all_y Logical. Whether or not to plot all plots on the same y axis limits.
#'@param y_log_axis Logical. Whether or not to put y axis on log scale
#'@param n_bins Number of bins for histograms. Default is 30.
#'@param n_cols Number of columns for faceted histograms. If NULL (default) will automatically choose n_cols for facetting.
#'@param text_size Size of text.
#'@param your_title Character. The title for the plot.
#'
#'@return Histogram of chosen statistics
#'
#'@importFrom magrittr %>%
#'@examples
#'stat_hist(your_SE = wu_subset[,1], data_choice = "assay stats",
#'          assay_choice = "counts")
#'@export
stat_hist <- function(your_SE,
                      data_choice = "assay stats",
                      assay_choice = "counts",
                      metadata_stat = NULL,
                      group_meta_by = NULL,
                      scale_all_y = FALSE,
                      y_log_axis = FALSE,
                      text_size = 12,
                      n_bins = 30,
                      n_cols = NULL,
                      your_title = NULL){


  if (data_choice == "assay stats"){

    # Error handling
    if (assay_choice %in% names(assays(your_SE)) == FALSE){
      stop("chosen assay is not in your_SE.")
    }

    # Load data
    your_data <- SummarizedExperiment::assays(your_SE)[[assay_choice]]

    # Make plots
    plot_list <- lapply(seq_len(ncol(your_data)), function(i){
      if (is.null(your_title)){
        your_title <- colnames(your_data)[i]
      }

      g <-ggplot2::ggplot(your_data, ggplot2::aes(x = your_data[,i]))+
        ggplot2::geom_histogram(bins = n_bins, color = "white", fill = "dodgerblue2")+
        cowplot::theme_cowplot()+
        ggplot2::labs(x = paste0("barcode ",assay_choice))+
        ggplot2::ggtitle(your_title)+
        ggplot2::theme(text = ggplot2::element_text(size = text_size), plot.title = ggplot2::element_text(face = "plain"))

    })

    if(scale_all_y){
      message((unlist(lapply(plot_list, function(x){ggplot2::ggplot_build(x)$layout$panel_params[[1]]$y.range[2]}))))
      plot_max <- max(unlist(lapply(plot_list, function(x){ggplot2::ggplot_build(x)$layout$panel_params[[1]]$y.range[2]})))
      message(plot_max)
      plot_list <- lapply(plot_list, function(x){x + ggplot2::coord_cartesian(ylim = c(NA, plot_max))})
    }

    if(y_log_axis){
      plot_list <- lapply(plot_list, function(x){x + ggplot2::scale_y_continuous(trans = "log10")})
    }

    if(is.null(n_cols)){
      g <- cowplot::plot_grid(plotlist = plot_list)
    } else {
      g <- cowplot::plot_grid(plotlist = plot_list, ncol = n_cols)
    }

  } else if (data_choice == "metadata stats"){

    # Load metadata
    meta_data <- as.data.frame(SummarizedExperiment::colData(your_SE))

    # Error handling
    if (metadata_stat %in% colnames(meta_data) == FALSE){
      stop("metadata_stat is not a piece of colData.")
    }

    if (is.null(group_meta_by)){
      meta_data_for_plot <- meta_data[,metadata_stat,drop = FALSE] %>%
        tibble::rownames_to_column(var = "samplename") %>%
        dplyr::rename(metadata_col = all_of(metadata_stat))
      if(is.numeric(meta_data_for_plot$metadata_col)){
        meta_data_for_plot$metadata_col <- as.factor(meta_data_for_plot$metadata_col)
      }

      g <- ggplot2::ggplot(meta_data_for_plot, ggplot2::aes(x = .data$metadata_col))+
        ggplot2::geom_histogram(bins = n_bins, position = "identity", stat = "count", fill = "dodgerblue2")+
        cowplot::theme_cowplot()+
        ggplot2::labs(x = paste0("metadata: ", metadata_stat)) +
        ggplot2::guides(color = FALSE)+
        ggplot2::theme(text = ggplot2::element_text(size = text_size), axis.text.x = element_text(angle = 45))
      if(y_log_axis){
        g <- g + ggplot2::scale_y_continuous(trans = "log10")
      }

    } else {
      if (!(group_meta_by %in% colnames(meta_data))){
        stop("group_meta_by is not a piece of colData.")
      }

      if(metadata_stat == group_meta_by){
        stop("cannot have metadata_stat and group_meta_by be the same column in colData")
      }
      meta_data_for_plot <- meta_data[,c(metadata_stat, group_meta_by)] %>%
        tibble::rownames_to_column(var = "samplename") %>%
        dplyr::rename(grouping_var = all_of(group_meta_by), metadata_col = all_of(metadata_stat))

      if(is.numeric(meta_data_for_plot$metadata_col)){
        meta_data_for_plot$metadata_col <- as.factor(meta_data_for_plot$metadata_col)
      }

      g <- ggplot2::ggplot(meta_data_for_plot, ggplot2::aes(x = .data$metadata_col))+
        ggplot2::geom_histogram(bins = n_bins, stat = "count",fill = "dodgerblue2")+
        cowplot::theme_cowplot()+
        ggplot2::labs(x = paste0("metadata: ", metadata_stat)) +
        ggplot2::guides(color = FALSE)+
        ggplot2::theme(text = ggplot2::element_text(size = text_size), axis.text.x = element_text(angle = 45))+
        ggplot2::facet_wrap(~grouping_var)

      if(y_log_axis){
        g <- g + ggplot2::scale_y_continuous(trans = "log1p")
      }
    }

  } else {
    stop("data_choice must be one of 'assay stats' or 'metadata stats' .")
  }

  g


}

