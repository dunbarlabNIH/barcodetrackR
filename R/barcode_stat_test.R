#' Barcode Statistical Test
#'
#' Carries out a specific instance of statistical testing relevant to clonal tracking experiments. For longitudinal observations (of barcode abundances) in the provided SE object, use a Chi-squared or Fisher exact test whether each barcode proportion has changed between samples. \cr Each column in the provided SE will be "tested" against the reference sample. If the `stat_option` argument is set to its default of "subsequent" then each sample will be compared to the sample before it. If this argument is set to "reference" the reference sample column name must be provided and each column will be tested against that reference sample.
#'
#'@param your_SE A Summarized Experiment object containing clonal tracking data as created by the barcodetrackR `create_SE` function.
#'@param sample_size A numeric vector providing the sample size of each column of the SummarizedExperiment passed to the function. This sample size describes the samples that the barcoding data is meant to approximate, for example the number of cells barcodes were extracted from.
#'@param stat_test The statistical test to use on the constructed contingency table for each barcode. Options are "chi-squared" and "fisher." \cr For information, see [chisq.test](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/chisq.test) [fisher.test](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/fisher.test)
#'@param stat_option For "subsequent" statistical testing is performed on each column of data compared to the column before it. For "reference," all other columns of data are compared to a reference column specified in the `reference_sample` arguument.
#'@param reference_sample Provide the column name of the reference column if stat_option is set to "reference." Defaults to the first column in the SummarizedExperiment.
#'@param p_adjust Character, default = "none". To correct p-values for muiltiple comparisons, set to any of the p value adjustment methods in the p.adjust function in R stats, which includes "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", and "fdr".
#'@param bc_threshold Clones must be above this proportion in at least one sample to be included in statistical testing. Default is 0. Use this to ignore low-abundance clones which are more likely to be noise or artifact.
#'@return Returns a list of 3 dataframes containing the following information for each observation (or barcode) which passed the provided bc_threshold: \cr [["FC"]], Fold Change of barcode abundance for each sample relative to the previous sample or to the specified reference sample. Please note that for maximal user control over results, the FC dataframe will contain 0 for barcodes where the test sample has an abundance of 0, Inf for barcodes where the reference sample had an abundance of 0 and NaN for a barcode where both the test and reference sample have an abundance of 0; \cr [["log_FC"]], same as previous but the log Fold Change. Please note that again for maximal user control, the log_FC dataframe will contain NaN values when the FC was Nan, -Inf values when the FC was 0, and Inf values when the FC was Inf; \cr [["p_val"]], the p-value returned from either the Chi-squared or Fisher exact test indicating whether each barcode changed in proportion between the test sample and the reference sample. Please note that the p value will be NaN if both abundances are 0, otherwise a p-value will be assigned. \cr Also, note that one column of each resulting dataframe will contain all NAs - in the case where the `stat_option` argument is set to "subsequent" then this will be the first sample since there is no subsequent sample to compare to. In the case where the `stat_option` argument is set to "reference" then the reference sample will contain NAs.
#'
#'@importFrom rlang %||%
#'@importFrom stats chisq.test
#'
#'@export
#'
#'@examples
#'barcode_stat_test(your_SE = wu_subset[,1:4], sample_size = rep(5000,4),
#'                  stat_test = "chi-squared", stat_option = "subsequent",
#'                  bc_threshold = 0.0001)
#'
barcode_stat_test <- function(your_SE,
                              sample_size,
                              stat_test = "chi-squared",
                              stat_option = "subsequent",
                              reference_sample = NULL,
                              p_adjust = "none",
                              bc_threshold = 0) {

  # Apply bc_threshold
  bc_passing_threshold <- apply(SummarizedExperiment::assays(your_SE)$proportions, 1, function(x){any(x>bc_threshold, na.rm = TRUE)})
  your_SE <- your_SE[bc_passing_threshold,]

  #error checking
  if (stat_test != "chi-squared" & stat_test != "fisher"){
    stop("stat_test must be either 'chi-squared' or 'fisher'.")
  }

  # Initialize dataframes of the right dimensions
  FC_df <- SummarizedExperiment::assays(your_SE)$proportions
  log_FC_df <- SummarizedExperiment::assays(your_SE)$proportions
  p_val_df <- SummarizedExperiment::assays(your_SE)$proportions

  if (stat_option == "subsequent"){
    stat_ref_index <- 1
    stat_test_index <- 2:length(colnames(your_SE))

    # Fill the first sample with default values
    FC_df[,1] <- rep(NA, times = nrow(your_SE))
    log_FC_df[,1] <- rep(NA, times = nrow(your_SE))
    p_val_df[,1] <- rep(NA, times = nrow(your_SE))

    # loop through each sample
    for (i in stat_test_index){

      # Calculate FC and log_FC
      FC_df[,i] <- SummarizedExperiment::assays(your_SE)$proportions[,i] / SummarizedExperiment::assays(your_SE)$proportions[,i-1]
      log_FC_df[,i] <- log2(FC_df[,i])

      # Perform statistical test
      if (stat_test == "chi-squared"){
        p_val_df[,i] <- sapply(1:nrow(your_SE), function(z) chisq.test(x = rbind(c(SummarizedExperiment::assays(your_SE)$proportions[z,i]*sample_size[i],
                                                                                SummarizedExperiment::assays(your_SE)$proportions[z,i-1]*sample_size[i-1]),
                                                                              c(sample_size[i] - SummarizedExperiment::assays(your_SE)$proportions[z,i]*sample_size[i],
                                                                                sample_size[i-1] - SummarizedExperiment::assays(your_SE)$proportions[z,i-1]*sample_size[i-1])))$p.val)
      } else if (stat_test == "fisher") {
        p_val_df[,i] <- sapply(1:nrow(your_SE), function(z) fisher.test(x = rbind(c(SummarizedExperiment::assays(your_SE)$proportions[z,i]*sample_size[i],
                                                                                 SummarizedExperiment::assays(your_SE)$proportions[z,i-1]*sample_size[i-1]),
                                                                               c(sample_size[i] - SummarizedExperiment::assays(your_SE)$proportions[z,i]*sample_size[i],
                                                                                 sample_size[i-1] - SummarizedExperiment::assays(your_SE)$proportions[z,i-1]*sample_size[i-1])))$p.val)
      }
    }

  } else if (stat_option == "reference"){
    stat_ref_choice <- reference_sample %||% colnames(your_SE)[1]
    # Error checking
    if (stat_ref_choice %in% colnames(your_SE) == FALSE){
      stop("reference_sample must be a column name in your_SE")
    }
    stat_test_index <- 1:length(colnames(your_SE))
    stat_ref_index <- which(colnames(your_SE) %in% stat_ref_choice)
    stat_test_index <- stat_test_index[!stat_test_index %in% stat_ref_index] # Remove reference column

    # Fill the reference sample with default values
    FC_df[,stat_ref_index] <- rep(NA, times = nrow(your_SE))
    log_FC_df[,stat_ref_index] <- rep(NA, times = nrow(your_SE))
    p_val_df[,stat_ref_index] <- rep(NA, times = nrow(your_SE))

    for (i in stat_test_index){

      # Calculate FC and log_FC
      FC_df[,i] <- SummarizedExperiment::assays(your_SE)$proportions[,i] / SummarizedExperiment::assays(your_SE)$proportions[,i-1]
      log_FC_df[,i] <- log2(FC_df[,i])

      # Perform statistical test compared to reference
      if (stat_test == "chi-squared"){
        p_mat[,i] <- sapply(1:nrow(your_SE), function(z) chisq.test(x = rbind(c(SummarizedExperiment::assays(your_SE)$proportions[z,i]*sample_size[i],
                                                                                SummarizedExperiment::assays(your_SE)$proportions[z,stat_ref_index]*sample_size[stat_ref_index]),
                                                                              c(sample_size[i] - SummarizedExperiment::assays(your_SE)$proportions[z,i]*sample_size[i],
                                                                                sample_size[stat_ref_index] - SummarizedExperiment::assays(your_SE)$proportions[z,stat_ref_index]*sample_size[stat_ref_index])))$p.val)
      } else if (stat_test == "fisher") {
        p_mat[,i] <- sapply(1:nrow(your_SE), function(z) fisher.test(x = rbind(c(SummarizedExperiment::assays(your_SE)$proportions[z,i]*sample_size[i],
                                                                                 SummarizedExperiment::assays(your_SE)$proportions[z,stat_ref_index]*sample_size[stat_ref_index]),
                                                                               c(sample_size[i] - SummarizedExperiment::assays(your_SE)$proportions[z,i]*sample_size[i],
                                                                                 sample_size[stat_ref_index] - SummarizedExperiment::assays(your_SE)$proportions[z,stat_ref_index]*sample_size[stat_ref_index])))$p.val)
      }
    }
  }

  # Correct p values for multiple comparisons
  p_val_df_adj <- as.data.frame(apply(p_val_df, 2, function(x) stats::p.adjust(x, method = p_adjust)))

  # Compile results into list
  output_list <- list(FC = FC_df,
                      log_FC = log_FC_df,
                      p_val = p_val_df_adj)

  return(output_list)
}








