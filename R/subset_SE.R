#' @title subset_SE
#'
#' @description Subsets an existing SummarizedExperiment object.
#'
#' @param your_SE A SummarizedExperiment object.
#' @param ... Arguments passed to subset_SE in the form of `X = keys` where `X` is a column from SE's colData and `keys` are entries in the colData to subset.
#' @return Returns a subsetted SummarizedExperiment object.
#' @examples
#' data(wu_subset)
#' wu_B.5month <- subset_SE(wu_subset, celltype = "B", timepoint = "6.0")
#' @export
subset_SE <- function(your_SE, ...) {
    arguments <- list(...)
    if (length(arguments) > 0) {
        for (i in seq_along(arguments)) {
            subset_name <- names(arguments)[[i]]
            subset_vars <- arguments[[i]]
            your_SE <- your_SE[, your_SE[[subset_name]] %in% subset_vars]
        }
    }
    return(your_SE)
}
