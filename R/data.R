#' Small subset of Wu barcoding dataset
#'
#' A SummarizedExperiment object containing a subset of the Wu barcoding dataset.
#' It includes peripheral blood T, B, Gr, NK_56, and NK-16 samples from the first
#' 4 times points of macaque ZJ31.
#'
#'
#' @format A SummarizedExperiment object with 215 features rows and 20 samples:
#' \describe{
#'   \item{assays}{includes the counts, percentages, ranks, normalized, and logs assays}
#'   \item{colData}{includes the accompanying metadata for the samples}
#'   \item{metadata}{includes the scale_factor used and the log_base used in the log assay}
#'   ...
#' }
#' @source
#' system.file("sample_data/WuC_etal/monkey_ZJ31.txt", package = "barcodetrackR") %>% read.delim(row.names = 1) -> wu_dataframe
#' system.file("sample_data/WuC_etal/monkey_ZJ31_metadata.txt", package = "barcodetrackR") %>% read.delim() -> wu_metadata
#' wu_SE <- create_SE(your_data = wu_dataframe, meta_data = wu_metadata, threshold = 0.005)
#' wu_subset <- wu_SE[,1:20]
#' \url{http://dx.doi.org/10.1126/sciimmunol.aat9781}
"wu_subset"
