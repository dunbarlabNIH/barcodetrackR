#' Rank Abundance Plot
#'
#' Create a rank abundance plot of the barcodes in the chosen samples provided in `your_SE`.
#' Use this function to visualize the distribution of barcode abundances within sample(s).
#' Note: If comparing the visualization to the statistical testing results from `rank_abundance_stat_test` function in barcodetrackR, please set the `scale_rank` to TRUE. The K-S test is agnostic to number of samples so it is directly comparable to the visualization produced when the barcode ranks are scaled between 0 and 1.
#'
#' @param your_SE Summarized Experiment object containing clonal tracking data as created by the barcodetrackR `create_SE` function.
#' @param scale_rank Logical. Whether or not to scale all ranks from 0 to 1 or keep barcode ranks as their actual integer values. When `scale_rank` is set to FALSE, all samples will not necessarily have the same x maximum.
#' @param point_size Numeric. Size of the points for the plot.
#' @param your_title Character. Specify a title for the plot.
#' @param text_size Numeric. Size of text in plot.
#' @param plot_labels Vector of labels for each sample. If not specified, the colnames(your_SE) will be used.
#' @param return_table Logical. If set to TRUE, rather than a plot, the function will return a dataframe containing for each sample, each barcode in rank order with its abundance in that sample, its scaled rank (0 to 1), and the cumulative sum of abundance for all barcodes with rank <= the rank of that barcode.
#'
#' @return Displays a rank-abundance plot (made by ggplot2) of the samples chosen. \cr Each point represents a single barcode with the x-value describing its rank in abundance with 1 being the most abundant barcode \cr The y-value representing the cumulative abundance of all barcodes with rank less than or equal to the x-axis value. \cr If the return_table is set to TRUE, instead of a plot, a datframe with the rank abundance data will be returned.
#'
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' rank_abundance_plot(your_SE = wu_subset[, 1:4], point_size = 2)
rank_abundance_plot <- function(your_SE,
    scale_rank = FALSE,
    point_size = 3,
    your_title = NULL,
    text_size = 12,
    plot_labels = NULL,
    return_table = FALSE) {

    # get labels
    plot_labels <- plot_labels %||% colnames(your_SE)
    if (length(plot_labels) != ncol(your_SE)) {
        stop("plot_labels must be same length as number of columns being plotted")
    }
    colnames(your_SE) <- plot_labels

    your_data <- SummarizedExperiment::assays(your_SE)[["proportions"]]

    lapply(seq_len(ncol(your_data)), function(i) {
        tibble::tibble(sample_name = colnames(your_data)[i], percentage = your_data[, i]) %>%
            dplyr::filter(.data$percentage > 0) %>%
            dplyr::arrange(desc(.data$percentage)) %>%
            dplyr::mutate(cumulative_sum = cumsum(.data$percentage), rank = dplyr::row_number()) %>%
            dplyr::mutate(scaled_rank = dplyr::percent_rank(-.data$percentage))
    }) %>%
        do.call(rbind, .) %>%
        dplyr::mutate(sample_name = factor(.data$sample_name, levels = colnames(your_data))) -> plotting_data

    scale_rank_choice <- ifelse(scale_rank, "scaled_rank", "rank")

    if (return_table) {
        return(plotting_data)
    }

    ggplot2::ggplot(plotting_data, ggplot2::aes_string(x = scale_rank_choice, y = "cumulative_sum", group = "sample_name", color = "sample_name")) +
        ggplot2::geom_point(size = point_size) +
        ggplot2::scale_color_discrete(labels = plot_labels) +
        ggplot2::theme_classic() +
        ggplot2::ggtitle(your_title) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(), text = ggplot2::element_text(size = text_size))
}
