#' Rank abundance plot
#'
#' Rank abundance Plot of the Barcodes for each sample chosen. Uses the cumulative sum so that all
#' samples eventually reach 100%.
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param point_size Numeric. Size of the points.
#'@param your_title The title for the plot.
#'@param scale_rank Whether or not to scale all ranks to 1 to 100 or keep numerical integer ranks.
#'@return Displays a rank-abundance plot (made by ggplot2) of the samples chosen.
#'
#'@importFrom magrittr %>%
#'
#'@export
#'
#'@examples
#'rank_abundance_plot(your_data = wu_SE, point_size = 4)
rank_abundance_plot = function(your_SE,
                               point_size = 3,
                               your_title = NULL,
                               scale_rank = FALSE) {

  your_data <- SummarizedExperiment::assays(your_SE)[["percentages"]]

  lapply(1:ncol(your_data), function(i){
    tibble::tibble(sample_name = colnames(your_data)[i], percentage = your_data[,i]) %>%
      dplyr::filter(percentage > 0) %>%
      dplyr::arrange(desc(percentage)) %>%
      dplyr::mutate(cumulative_sum = cumsum(percentage), rank = dplyr::row_number()) %>%
      dplyr::mutate(scaled_rank = dplyr::percent_rank(-percentage))
  }) %>% do.call(rbind, .) %>% dplyr::mutate(sample_name = factor(sample_name, levels = colnames(your_data))) -> plotting_data

  scale_rank_choice <- ifelse(scale_rank, "scaled_rank", "rank")

  ggplot2::ggplot(plotting_data, ggplot2::aes_string(x = scale_rank_choice, y = "cumulative_sum", group = "sample_name", color = "sample_name"))+
    ggplot2::geom_point(size = point_size)+
    ggplot2::theme_bw()+
    ggplot2::ggtitle(your_title)+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), panel.grid = ggplot2::element_blank())
}
