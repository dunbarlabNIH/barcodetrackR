context("Shared Clonality Functions")

# Load sample data
test_data <- read.delim(system.file("sample_data/app_sample_data/sample_data_ZJ31.txt", package = "barcodetrackR"), row.names = 1)
test_metadata <- read.delim(system.file("sample_data/app_sample_data/sample_metadata_ZJ31.txt", package = "barcodetrackR"))
my_SE <- barcodetrackR::create_SE(your_data = test_data,
                                  meta_data  = test_metadata)

test_that("scatter_plot works", {
  testthat::expect_type(barcodetrackR::scatter_plot(my_SE[,1:2]), "list")
})
#> Test passed 🥳

test_that("cor_plot works", {
  testthat::expect_type(barcodetrackR::cor_plot(my_SE[,1:3]), "list")
  testthat::expect_s3_class(barcodetrackR::cor_plot(my_SE[,1:3], return_table = TRUE), "data.frame")
})
#> Test passed 🥳

test_that("chord_diagram works", {
  testthat::expect_type(barcodetrackR::chord_diagram(my_SE[,1:3]), "list")
  testthat::expect_s3_class(barcodetrackR::chord_diagram(my_SE[,1:3], return_table = TRUE), "data.frame")
})
#> Test passed 🥳

test_that("autocor_plot works", {
  testthat::expect_type(barcodetrackR::autocor_plot(my_SE,
                                                    filter_by = "celltype",
                                                    filter_selection = "T",
                                                    plot_over = "months"), "list")
  testthat::expect_s3_class(barcodetrackR::autocor_plot(my_SE,
                                                        filter_by = "celltype",
                                                        filter_selection = "T",
                                                        plot_over = "months",
                                                        return_table = TRUE), "data.frame")
})
#> Test passed 🥳