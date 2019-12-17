# diego espinoza
# 12/16/2019

# Purpose: to load barcode data and associated metadata into SummarizedExperiment
# to call our different visualization functions on data stored in this object

require(magrittr)

# Load barcode and meta data
barcode.file <- read.delim("~/Documents/dunbar_lab/projects/diego_barcodetrackR/barcodetrackR/data/ZG66_simple_data.txt", row.names = 1)
meta.file <- read.delim("~/Documents/dunbar_lab/projects/diego_barcodetrackR/barcodetrackR/data/ZG66_simple_metadata.txt") %>% dplyr::mutate(SAMPLENAME = FILENAME)

# Create summarized experiment
SE <- create_SE(your_data = barcode.file, meta_data = meta.file)

# Heat map of just the B cells
barcodetrackR::barcode_ggheatmap(your_data = SummarizedExperiment::assay(se)[,se$Cell_type == "B"],
                                 names = SummarizedExperiment::colData(se[,se$Cell_type == "B"])$Timepoint,
                                 n_clones = 10, your_title = "ZG66 B cells", label_size = 12)

# Ridge plot of T cells vs B cells
source('ridge_plot.R')
source('ridge_plot_weighted.R')
ridge_plot(your_data = SummarizedExperiment::assay(se))
ridge_plot_weighted(your_data = SummarizedExperiment::assay(se))

# Correlation plot
barcodetrackR::cor_plot(your_data = SummarizedExperiment::assay(se),
                        names = SummarizedExperiment::colData(se)$FILENAME,
                        your_title = "Correlation plot",
                        labelsizes = 1, plottype = "ellipse")

# Calculate diversity
source('clonal_diversity.R')
my.div <- clonal_diversity(your_data = SummarizedExperiment::assay(se), divindex = "shannon")

# Diversity plot
source('diversity_plot.R')
# Plot diversity
diversity_plot(your_SE = se, index_type = "shannon", measure = "diversity", plot_by = "Timepoint", group_by = "Cell_type")
# Plot Shannon count
diversity_plot(your_SE = se, index_type = "shannon", measure = "count", plot_by = "Timepoint", group_by = "Cell_type")
# Plot Shannon evenness
diversity_plot(your_SE = se, index_type = "shannon", measure = "evenness", plot_by = "Timepoint", group_by = "Cell_type")

