context("Clonal Diversity Functions")

# Load sample data
test_data <- read.delim(system.file("sample_data/app_sample_data/sample_data_ZJ31.txt", package = "barcodetrackR"), row.names = 1)
test_metadata <- read.delim(system.file("sample_data/app_sample_data/sample_metadata_ZJ31.txt", package = "barcodetrackR"))
my_SE <- barcodetrackR::create_SE(your_data = test_data,
                                  meta_data  = test_metadata)

test_that("rank_abundance works", {
  testthat::expect_type(barcodetrackR::rank_abundance_plot(my_SE[,1:5]), "list")
  testthat::expect_s3_class(barcodetrackR::rank_abundance_plot(my_SE[,1:5], return_table = TRUE), "data.frame")
})
#> Test passed 🥳

test_that("clonal_diversity works", {
  testthat::expect_type(barcodetrackR::clonal_diversity(my_SE,
                                                        plot_over = "months",
                                                        group_by = "celltype"), "list")
  testthat::expect_s3_class(barcodetrackR::clonal_diversity(my_SE,
                                                            plot_over = "months",
                                                            group_by = "celltype",
                                                            return_table = TRUE), "data.frame")
})
#> Test passed 🥳

test_that("clonal_count works", {
  testthat::expect_type(barcodetrackR::clonal_count(my_SE,
                                                    plot_over = "months",
                                                    group_by = "celltype"), "list")
  testthat::expect_s3_class(barcodetrackR::clonal_count(my_SE,
                                                        plot_over = "months",
                                                        group_by = "celltype",
                                                        return_table = TRUE), "data.frame")
})
#> Test passed 🥳

test_that("mds_plot works", {
  testthat::expect_type(barcodetrackR::mds_plot(my_SE[,1:5]), "list")
  testthat::expect_s3_class(barcodetrackR::mds_plot(my_SE[,1:5], return_table = TRUE), "data.frame")
})
#> Test passed 🥳