# Testing barcodetrackR code

system.file("sample_data/WuC_etal/monkey_ZJ31.txt", package = "barcodetrackR") %>%
  read.delim(row.names = 1) -> wu_dataframe

system.file("sample_data/WuC_etal/monkey_ZJ31_metadata.txt", package = "barcodetrackR") %>%
  read.delim() -> wu_metadata

wu_SE <- barcodetrackR::create_SE(your_data = wu_dataframe,
                   meta_data = wu_metadata,
                   threshold = 0.005)

# Test barcode_statistics
my_stats <- barcode_statistics(your_SE = wu_SE[,1:20],
                               sample_size = rep(5000, 20))
