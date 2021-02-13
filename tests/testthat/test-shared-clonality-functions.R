context("Shared Clonality Functions")

test_that("scatter_plot works", {
  testthat::expect_type(barcodetrackR::scatter_plot(wu_subset[,1:2]), "list")
})
#> Test passed 🥳

test_that("cor_plot works", {
  testthat::expect_type(barcodetrackR::cor_plot(wu_subset[,1:3]), "list")
  testthat::expect_s3_class(barcodetrackR::cor_plot(wu_subset[,1:3], return_table = TRUE), "data.frame")
})
#> Test passed 🥳

test_that("chord_diagram works", {
  testthat::expect_type(barcodetrackR::chord_diagram(wu_subset[,1:3]), "list")
  testthat::expect_s3_class(barcodetrackR::chord_diagram(wu_subset[,1:3], return_table = TRUE), "data.frame")
})
#> Test passed 🥳

test_that("autocor_plot works", {
  testthat::expect_type(barcodetrackR::autocor_plot(wu_subset,
                                                    filter_by = "celltype",
                                                    filter_selection = "T",
                                                    plot_over = "months"), "list")
  testthat::expect_s3_class(barcodetrackR::autocor_plot(wu_subset,
                                                        filter_by = "celltype",
                                                        filter_selection = "T",
                                                        plot_over = "months",
                                                        return_table = TRUE), "data.frame")
})
#> Test passed 🥳