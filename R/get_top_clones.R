#' @title get_top_clones (helper function)
#'
#' @description Retrieves the sequence(s) (row-identifier(s)) of the top "n_clones" within the specified sample from a SummarizedExperiment object.
#'
#' @param your_SE A summarized experiment.
#' @param SAMPLENAME_choice Name of the SAMPLENAME identifier within your_SE from which to retrieve the top clones from.
#' @param n_clones Numeric. Number of top clones from the specified sample that should be retrieved.
#'
#' @return The row indices for the top n_clones in the dataset, using the 'ranks' assay.
#' @examples
#' get_top_clones(wu_subset, "ZJ31_6m_T", n_clones = 10)
#' @export

get_top_clones <- function(your_SE, SAMPLENAME_choice, n_clones = 10) {
    sample_choice <- SummarizedExperiment::assays(your_SE)[["ranks"]][, SAMPLENAME_choice, drop = FALSE]
    indices <- rownames(sample_choice)[which(sample_choice <= n_clones)]
    return(indices)
}
