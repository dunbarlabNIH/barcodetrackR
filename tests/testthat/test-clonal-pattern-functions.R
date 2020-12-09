context("Clonal Pattern Functions")

# Load sample data
test_data <- read.delim(system.file("sample_data/app_sample_data/sample_data_ZJ31.txt", package = "barcodetrackR"), row.names = 1)
test_metadata <- read.delim(system.file("sample_data/app_sample_data/sample_metadata_ZJ31.txt", package = "barcodetrackR"))
my_SE <- barcodetrackR::create_SE(your_data = test_data,
                                  meta_data  = test_metadata)

test_that("barcode_ggheatmap works", {
  testthat::expect_type(barcodetrackR::barcode_ggheatmap(my_SE), "list")
  testthat::expect_s3_class(barcodetrackR::barcode_ggheatmap(my_SE, return_table = TRUE), "data.frame")
})
#> Test passed 🥳

test_that("barcode_ggheatmap_stat works", {
  testthat::expect_type(barcodetrackR::barcode_ggheatmap_stat(my_SE[1:10000,1:3],
                                                              sample_size = rep(1000, ncol(my_SE))), "list")
  testthat::expect_s3_class(barcodetrackR::barcode_ggheatmap_stat(my_SE[1:10000,1:3],
                                                                  sample_size = rep(1000, ncol(my_SE)),
                                                                  return_table = TRUE), "data.frame")
})
#> Test passed 🥳

test_that("barcode_binary_heatmap works", {
  testthat::expect_type(barcodetrackR::barcode_binary_heatmap(my_SE), "list")
  testthat::expect_s3_class(barcodetrackR::barcode_binary_heatmap(my_SE, return_table = TRUE), "data.frame")
})
#> Test passed 🥳

test_that("clonal_contribution works", {
  testthat::expect_type(barcodetrackR::clonal_contribution(my_SE,
                                                           SAMPLENAME_choice = "ZJ31_20m_T",
                                                           filter_by = "celltype",
                                                           filter_selection = "T",
                                                           plot_over = "months"), "list")
  testthat::expect_s3_class(barcodetrackR::clonal_contribution(my_SE,
                                                               SAMPLENAME_choice = "ZJ31_20m_T",
                                                               filter_by = "celltype",
                                                               filter_selection = "T",
                                                               plot_over = "months",
                                                               return_table = TRUE), "data.frame")
})
#> Test passed 🥳