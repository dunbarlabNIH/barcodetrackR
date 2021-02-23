#' MDS Plot
#'
#' Calculates a simmilarity/dissimlarity index or metrix for each sample-sample pair and reduces the resulting dist matrix into two dimensions
#'
#'@param your_SE Summarized Experiment object containing clonal tracking data as created by the barcodetrackR `create_SE` function.
#'@param group_by Column of metadata to color samples by. Can also specify "kmeans_cluster" if kmeans_cluster argument is set to TRUE, and then the grouping variables will be the clusterinng result.
#'@param method_dist Dissimilarity index from vegan. One of "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", or "cao".
#'@param assay The assay to calculate the index on
#'@param your_title Character. The title for the plot.
#'@param point_size Numeric. The size of the points.
#'@param text_size Numeric. Size of text in plot.
#'@param return_table Logical. If set to true, the function will return a dataframe containing each samples reduced measure of dissimilarity coordinates.
#'@param kmeans_cluster Logical. If set to true, each sample will be assigned a cluster computed by kmeans on the chosen assay.
#'@param k.param Numeric. If kmeans_cluster is TRUE, provide the number of kmeans clusters to identify.
#'@param draw_ellipses Logical. If kmeans_cluster is TRUE, draw ellipses around the different kmeans clusters.
#'
#'@return Plots dissimilarity indices between samples in your_SE. Or if return table is set to TRUE, returns a dataframe of each sample's reduced measures of dissimilarity coordinates.
#'
#'@importFrom rlang %||%
#'@importFrom magrittr %>%
#'
#'@export
#'
#'@examples
#'mds_plot(your_SE = wu_subset, method_dist = "bray",group_by = "celltype")
#"

mds_plot = function(your_SE,
                    group_by = "SAMPLENAME",
                    method_dist ="bray",
                    assay = "proportions",
                    your_title = NULL,
                    point_size = 3,
                    text_size = 12,
                    return_table = FALSE,
                    kmeans_cluster = FALSE,
                    k.param = 3,
                    draw_ellipses = FALSE){

  SummarizedExperiment::colData(your_SE) %>%
    tibble::as_tibble() %>%
    dplyr::mutate_if(is.factor, as.character) -> your_colData

  #extracts chosen assay from your_SE
  t(SummarizedExperiment::assays(your_SE)[[assay]]) %>%
    vegan::vegdist(method = method_dist) %>%
    stats::cmdscale() %>%
    magrittr::set_colnames(c("MDS_1", "MDS_2")) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    tibble::rownames_to_column(var = "SAMPLENAME") %>%
    dplyr::left_join(your_colData, by = "SAMPLENAME") -> plotting_data

  # if(is.null(your_title)){
  #   your_title <- paste0("PCoA on pairwise ", method_dist, " values")
  # }

  if (kmeans_cluster){
    plotting_data$kmeans_cluster <- as.factor(stats::kmeans(t(SummarizedExperiment::assays(your_SE)[[assay]]), k.param)$cluster)
    # plotting_data$kmeans_cluster <- as.factor(kmeans(plotting_data[,c("MDS_1","MDS_2")], k.param)$cluster) # if we wanted to cluster in mds space.
    if (draw_ellipses){
      if (min(table(plotting_data$kmeans_cluster)) < 4){
        stop("Please choose a lower k.param value. Ellipses cannot be drawn if less than 4 observations are in one k mean cluster.")
      }
    }
  }

  if (return_table){
    return(plotting_data)
  }

  p <- ggplot2::ggplot(plotting_data, ggplot2::aes_string(x = "MDS_1", y = "MDS_2", color = group_by)) +
    ggplot2::geom_point(size = point_size)+
    ggplot2::ggtitle(your_title)+
    ggplot2::theme_classic()+
    ggplot2::theme(text = ggplot2::element_text(size=text_size))

  if (kmeans_cluster & draw_ellipses){
    p <- p + ggplot2::stat_ellipse(aes(x = .data$MDS_1, y = .data$MDS_2, group = kmeans_cluster), color = "black", linetype = 2)
  }

  p
}






