context("Load data")

# Load sample data
test_data <- read.delim(system.file("sample_data/app_sample_data/sample_data_ZJ31.txt", package = "barcodetrackR"), row.names = 1)
test_metadata <- read.delim(system.file("sample_data/app_sample_data/sample_metadata_ZJ31.txt", package = "barcodetrackR"))

test_that("Creating SE works", {
  expect_s4_class(barcodetrackR::create_SE(your_data = test_data,
                                           meta_data  = test_metadata), "SummarizedExperiment")
})
#> Test passed 🥳

# Set threshold
test_threshold <- 0.0000001
test_that("Thresholding data works", {
  expect_s4_class(barcodetrackR::create_SE(your_data = test_data,
                                           meta_data  = test_metadata,
                                           threshold = test_threshold), "SummarizedExperiment")
})
#> Test passed 🥳

