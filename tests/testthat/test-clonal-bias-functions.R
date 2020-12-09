context("Clonal Bias Functions")

# Load sample data
test_data <- read.delim(system.file("sample_data/app_sample_data/sample_data_ZJ31.txt", package = "barcodetrackR"), row.names = 1)
test_metadata <- read.delim(system.file("sample_data/app_sample_data/sample_metadata_ZJ31.txt", package = "barcodetrackR"))
my_SE <- barcodetrackR::create_SE(your_data = test_data,
                                  meta_data  = test_metadata)

test_that("bias_histogram works", {
  testthat::expect_type(barcodetrackR::bias_histogram(my_SE,
                                                      split_bias_on = "celltype",
                                                      bias_1 = "B",
                                                      bias_2 = "Gr",
                                                      split_bias_over = "months"), "list")
  testthat::expect_type(barcodetrackR::bias_histogram(my_SE,
                                                          split_bias_on = "celltype",
                                                          bias_1 = "B",
                                                          bias_2 = "Gr",
                                                          split_bias_over = "months",
                                                          return_table = TRUE), "list")
})
#> Test passed 🥳

test_that("bias_ridge_plot works", {
  testthat::expect_type(barcodetrackR::bias_ridge_plot(my_SE,
                                                       split_bias_on = "celltype",
                                                       bias_1 = "B",
                                                       bias_2 = "Gr",
                                                       split_bias_over = "months"), "list")
  testthat::expect_s3_class(barcodetrackR::bias_ridge_plot(my_SE,
                                                           split_bias_on = "celltype",
                                                           bias_1 = "B",
                                                           bias_2 = "Gr",
                                                           split_bias_over = "months",
                                                           return_table = TRUE), "data.frame")
})
#> Test passed 🥳

test_that("bias_lineplot works", {
  testthat::expect_type(barcodetrackR::bias_lineplot(my_SE,
                                                     split_bias_on = "celltype",
                                                     bias_1 = "B",
                                                     bias_2 = "Gr",
                                                     split_bias_over = "months",
                                                     remove_unique = TRUE), "list")
  testthat::expect_s3_class(barcodetrackR::bias_lineplot(my_SE,
                                                         split_bias_on = "celltype",
                                                         bias_1 = "B",
                                                         bias_2 = "Gr",
                                                         split_bias_over = "months",
                                                         remove_unique = TRUE, 
                                                         return_table = TRUE), "data.frame")
})
#> Test passed 🥳