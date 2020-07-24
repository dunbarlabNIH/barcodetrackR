#' MDS Plot
#'
#' Calculates a simmilarity/dissimlarity index or metrix for each sample-sample pair and reduces the resulting dist matrix into two dimensions
#'
#'@param your_SE A Summarized Experiment object.
#'@param group_by Selection from colData used to label each point. Defaults to SAMPLENAME.
#'@param method_dist Dissimilarity index from vegan. One of "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", or "cao".
#'@param assay The assay to calculate the index on
#'@param your_title The title for the plot.
#'@param point_size The size of the points.
#'@param text_size Numeric. Size of text in plot.
#'@return Plots pairwise correlation plot for the samples in your_SE.
#'
#'@importFrom rlang %||%
#'@importFrom magrittr %>%
#'
#'@export
#'
#'@examples
#'mds_plot(your_SE = wu_SE, your_title = "MDS plot of Bray-Curtis dissimilarities", method_dist = "bray",group_by = "celltype")
#"

mds_plot = function(your_SE,
                    group_by = "SAMPLENAME",
                    method_dist ="bray",
                    assay = "percentages",
                    your_title = NULL,
                    point_size = 3,
                    text_size = 12) {

  SummarizedExperiment::colData(your_SE) %>%
    tibble::as_tibble() %>%
    dplyr::mutate_if(is.factor, as.character) -> your_colData

  #extracts chosen assay from your_SE
  t(SummarizedExperiment::assays(your_SE)[[assay]]) %>%
    vegan::vegdist(method = method_dist) %>%
    cmdscale() %>%
    magrittr::set_colnames(c("MDS_1", "MDS_2")) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    tibble::rownames_to_column(var = "SAMPLENAME") %>%
    dplyr::left_join(your_colData, by = "SAMPLENAME") -> plotting_data

  # if(is.null(your_title)){
  #   your_title <- paste0("PCoA on pairwise ", method_dist, " values")
  # }

  ggplot2::ggplot(plotting_data, ggplot2::aes_string(x = "MDS_1", y = "MDS_2", color = group_by)) +
    ggplot2::geom_point(size = point_size)+
    ggplot2::ggtitle(your_title)+
    ggplot2::theme_classic()+
    ggplot2::theme(text = ggplot2::element_text(size=text_size))

}






