#' Barcode Statistics
#'
#' Calculates statistical testing comparing each barcode between samples.
#'
#'@param your_SE A Summarized Experiment object.
#'@param sample_size A numeric vector providing the sample size of each column of the SummarizedExperiment passed to the function. This sample size describes the samples that the barcoding data is meant to approximate.
#'@param stat_test The statistical test to use on the constructed contingency table for each barcoe. Options are "chi-squared" and "fisher."
#'@param stat_option For "subsequent" statistical testing is performed on each column of data compared to the column before it. For "reference," all other columns of data are compared to a reference column.
#'@param reference_sample Provide the column name of the reference column if stat_option is set to "reference." Defaults to the first column in the SummarizedExperiment.
#'@param bc_threshold Clones must be above this proportion in at least one sample to be included in statistical testing.
#'@return Displays a heatmap in the current plot window.
#'
#'@importFrom rlang %||%
#'
#'@export
#'
#'@examples
#'barcode_statistics(your_SE = ZH33_SE, sample_size = rep(5000,ncol(ZH33_SE)), stat_test = "chi-squared", stat_option = "subsequent", bc_threshold = 0.0001)
#'
barcode_statistics <- function(your_SE,
                               sample_size,
                               stat_test = "chi-squared",
                               stat_option = "subsequent",
                               reference_sample = NULL,
                               bc_threshold = 0) {

  # Apply bc_threshold
  bc_passing_threshold <- apply(SummarizedExperiment::assays(your_SE)$percentages, 1, function(x){any(x>bc_threshold, na.rm = TRUE)})
  your_SE <- your_SE[bc_passing_threshold,]

  #error checking
  if (stat_test != "chi-squared" & stat_test != "fisher"){
    stop("stat_test must be either 'chi-squared' or 'fisher' for now.")
  }

  # Initialize dataframes of the right dimensions
  FC_df <- SummarizedExperiment::assays(your_SE)$percentages
  log_FC_df <- SummarizedExperiment::assays(your_SE)$percentages
  p_val_df <- SummarizedExperiment::assays(your_SE)$percentages

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
      FC_df[,i] <- SummarizedExperiment::assays(your_SE)$percentages[,i] / SummarizedExperiment::assays(your_SE)$percentages[,i-1]
      log_FC_df[,i] <- log2(FC_df[,i])

      # Perform statistical test
      if (stat_test == "chi-squared"){
        p_val_df[,i] <- sapply(1:nrow(your_SE), function(z) chisq.test(x = rbind(c(SummarizedExperiment::assays(your_SE)$percentages[z,i]*sample_size[i],
                                                                                SummarizedExperiment::assays(your_SE)$percentages[z,i-1]*sample_size[i-1]),
                                                                              c(sample_size[i] - SummarizedExperiment::assays(your_SE)$percentages[z,i]*sample_size[i],
                                                                                sample_size[i-1] - SummarizedExperiment::assays(your_SE)$percentages[z,i-1]*sample_size[i-1])))$p.val)
      } else if (stat_test == "fisher") {
        p_val_df[,i] <- sapply(1:nrow(your_SE), function(z) fisher.test(x = rbind(c(SummarizedExperiment::assays(your_SE)$percentages[z,i]*sample_size[i],
                                                                                 SummarizedExperiment::assays(your_SE)$percentages[z,i-1]*sample_size[i-1]),
                                                                               c(sample_size[i] - SummarizedExperiment::assays(your_SE)$percentages[z,i]*sample_size[i],
                                                                                 sample_size[i-1] - SummarizedExperiment::assays(your_SE)$percentages[z,i-1]*sample_size[i-1])))$p.val)
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
      FC_df[,i] <- SummarizedExperiment::assays(your_SE)$percentages[,i] / SummarizedExperiment::assays(your_SE)$percentages[,i-1]
      log_FC_df[,i] <- log2(FC_df[,i])

      # Perform statistical test compared to reference
      if (stat_test == "chi-squared"){
        p_mat[,i] <- sapply(1:nrow(your_SE), function(z) chisq.test(x = rbind(c(SummarizedExperiment::assays(your_SE)$percentages[z,i]*sample_size[i],
                                                                                SummarizedExperiment::assays(your_SE)$percentages[z,stat_ref_index]*sample_size[stat_ref_index]),
                                                                              c(sample_size[i] - SummarizedExperiment::assays(your_SE)$percentages[z,i]*sample_size[i],
                                                                                sample_size[stat_ref_index] - SummarizedExperiment::assays(your_SE)$percentages[z,stat_ref_index]*sample_size[stat_ref_index])))$p.val)
      } else if (stat_test == "fisher") {
        p_mat[,i] <- sapply(1:nrow(your_SE), function(z) fisher.test(x = rbind(c(SummarizedExperiment::assays(your_SE)$percentages[z,i]*sample_size[i],
                                                                                 SummarizedExperiment::assays(your_SE)$percentages[z,stat_ref_index]*sample_size[stat_ref_index]),
                                                                               c(sample_size[i] - SummarizedExperiment::assays(your_SE)$percentages[z,i]*sample_size[i],
                                                                                 sample_size[stat_ref_index] - SummarizedExperiment::assays(your_SE)$percentages[z,stat_ref_index]*sample_size[stat_ref_index])))$p.val)
      }
    }
  }

  # Compile results into list
  output_list <- list(FC = FC_df,
                      log_FC = log_FC_df,
                      p_val = p_val_df)

  return(output_list)
}








