#' Correlation Plot
#'
#' Gives the pairwise correlation between each sample-sample pair in the data frame.
#'
#'@param your_SE A Summarized Experiment object.
#'@param plot_labels Vector of x axis labels. Defaults to colnames(your_SE).
#'@param method_corr Character. One of "pearson", "spearman", or "kendall". Can also use "manhattan" to compute manhattan distance instead.
#'@param your_title The title for the plot.
#'@param grid Logical. Include a grid or not in the correlation plot
#'@param label_size The size of the column labels.
#'@param plot_type Character. One of "color", "circle", or "number".
#'@param no_negatives Logical. Whether to make negative correlations = 0.
#'@param return_table Logical. Whether or not to return table of p-values, confidence intervals, and R values instead of displaying a plot.
#'@param color_scale Character. Either "default" or an odd-numbered color scale where the lowest value will correspond to -1, the median value to 0, and the highest value to 1.
#'@param number_size Size of the text label when plot_type is "number".
#'@return Plots pairwise correlation plot for the samples in your_SE.
#'
#'@importFrom rlang %||%
#'@importFrom magrittr %>%
#'
#'@export
#'
#'@examples
#'cor_plot(your_SE = zh33, your_title = "Pearson correlation of all samples", plottype = "color")
#"

cor_plot = function(your_SE,
                    plot_labels = colnames(your_SE),
                    method_corr ="pearson",
                    your_title = "",
                    grid = TRUE,
                    label_size = 1,
                    plot_type = "color",
                    no_negatives = FALSE,
                    return_table = FALSE,
                    color_scale = "default",
                    number_size = 0.5) {

  #extracts percentages assay from your_SE
  plotting_data <- SummarizedExperiment::assays(your_SE)[["percentages"]]

  plotting_data_columns <- colnames(plotting_data)

  plotting_data_longer <- lapply(1:length(plotting_data_columns), function(i){
    lapply(1:length(plotting_data_columns), function(j){

      temp_df <- data.frame(plotting_data[[i]], plotting_data[[j]])
      temp_df <- temp_df[rowSums(temp_df) > 0, ]
      cortest_results <- cor.test(temp_df[[1]], temp_df[[2]], method = method_corr)
      result_df <- data.frame(sample_i = plot_labels[i],
                              sample_j = plot_labels[j],
                              correlation_value = cortest_results$estimate,
                              p_value = cortest_results$p.value,
                              ci_lo = cortest_results$conf.int[1],
                              ci_hi = cortest_results$conf.int[2])
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .) %>%
    dplyr::mutate(sample_i = factor(sample_i, levels = plot_labels), sample_j = factor(sample_j, levels = rev(plot_labels)))

  if(color_scale == "default"){
    color_scale <- c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")
  } else {
    if((length(color_scale) %% 2) != 1){
      stop("color_scale must be a vector of odd length")
    }
  }

  if(no_negatives){
    plotting_data_longer %>%
      dplyr::mutate(p_value = ifelse(correlation_value < 0, NA, p_value)) %>%
      dplyr::mutate(ci_lo = ifelse(correlation_value < 0, NA, ci_lo)) %>%
      dplyr::mutate(ci_hi = ifelse(correlation_value < 0, NA, ci_hi)) %>%
      dplyr::mutate(correlation_value = ifelse(correlation_value < 0, 0, correlation_value)) -> plotting_data_longer
    color_limits <- c(0, 1)
    floor_limit <- ceiling(length(color_scale)/2)
    color_scale <- color_scale[floor_limit:length(color_scale)]
  } else {
    color_limits <- c(-1,1)
  }

  if(return_table){
    return(plotting_data_longer)
  }

  gg_corplot <- ggplot2::ggplot(plotting_data_longer, ggplot2::aes(x = sample_i, y = sample_j)) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                   rect = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 0),
                   axis.title = ggplot2::element_blank())

  if(plot_type == "color"){
    gg_corplot <- gg_corplot +
      ggplot2::geom_tile(ggplot2::aes(fill = correlation_value), color = ifelse(grid, "black", NULL)) +
      ggplot2::scale_fill_gradientn(colours = color_scale, limits = color_limits, name = "correlation")
  } else if (plot_type == "circle"){
    gg_corplot <- gg_corplot +
      ggplot2::geom_tile(color = "black", fill = "white") +
      ggplot2::geom_point(ggplot2::aes(size = abs(correlation_value), fill = correlation_value), shape = 21)+
      ggplot2::scale_size_area("|correlation|", limits = c(0, 1)) +
      ggplot2::scale_fill_gradientn(colours = color_scale, limits = color_limits, name = "correlation")
  } else if (plot_type == "number") {
    gg_corplot <- gg_corplot +
      ggplot2::geom_tile(ggplot2::aes(fill = correlation_value), color = ifelse(grid, "black", NULL)) +
      ggplot2::scale_fill_gradientn(colours = color_scale, limits = color_limits, name = "correlation")+
      ggplot2::geom_text(ggplot2::aes(label = round(correlation_value, digits = 2)), color = "black", size = number_size)
  } else {
    stop("plot_type must be one of \"color\", \"circle\", or \"number\"")
  }

  return(gg_corplot)


}






