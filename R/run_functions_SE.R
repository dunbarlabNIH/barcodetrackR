
# Purpose: to load barcode data and associated metadata into SummarizedExperiment
# to call our different visualization functions on data stored in this object

require(tidyverse)

# Load barcode and meta data
barcode.file <- read.delim("/Users/mortlockrd/Desktop/GitHub/barcodetrackR/inst/sample_data/ZG66_simple_data.txt", row.names = 1)
meta.file <- read.delim("/Users/mortlockrd/Desktop/GitHub/barcodetrackR/inst/sample_data/ZG66_simple_metadata.txt") %>% dplyr::mutate(SAMPLENAME = FILENAME)

# Create summarized experiment
source('create_SE.R')
source('subset_SE.R')
source('threshold.R')
se <- create_SE(your_data = barcode.file, meta_data = meta.file, threshold = 5e-4)

# Heat map of just the B cells
barcodetrackR::barcode_ggheatmap(your_data = SummarizedExperiment::assay(se)[,se$Cell_type == "B"],
                                 names = SummarizedExperiment::colData(se[,se$Cell_type == "B"])$Timepoint,
                                 n_clones = 10, your_title = "ZG66 B cells", label_size = 12)

# Ridge plot of T cells vs B cells
source('ridge_plot.R')
# If plot_by is a factor, it will follow the order it appears in colData
ridge_plot(your_SE = se, cell_var = "Cell_type",cell_1 = "B", cell_2 = "T", plot_by = "Timepoint")
# If plot_by is numeric, it will sort from smallest to largest but still plot as categorical
ridge_plot(your_SE = se, cell_var = "Cell_type",cell_1 = "T", cell_2 = "B", plot_by = "Months", scale = 1.5)

# Weighted by overall contribution of each barcode
ridge_plot(your_SE = se, cell_var = "Cell_type",cell_1 = "B", cell_2 = "T", plot_by = "Timepoint", weighted = T)
ridge_plot(your_SE = se, cell_var = "Cell_type",cell_1 = "T", cell_2 = "B", plot_by = "Months", weighted = T,scale = 1.5)


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

# Clonal contribution
source('clonal_contribution.R')
clonal_contribution(your_SE = se, graph_type = "bar", filter_by = "Cell_type", filter_selection = "T", plot_by = "Timepoint", n_clones = 20, linesize = .4)
clonal_contribution(your_SE = se, graph_type = "line", filter_by = "Cell_type", filter_selection = "B", plot_by = "Timepoint", n_clones = 100, linesize = .1)

# Track top clones
source('track_top_clones.R')
track_top_clones(your_SE = se, filter_by = "Cell_type", filter_selection = "T", plot_by = "Timepoint", chosen_sample = "82hm", n_clones = 12, linesize = .4)
track_top_clones(your_SE = se, filter_by = "Cell_type", filter_selection = "B", plot_by = "Timepoint", chosen_sample = "82hm", n_clones = 25, linesize = .4)
track_top_clones(your_SE = se, filter_by = "Cell_type", filter_selection = "T", plot_by = "Timepoint", chosen_sample = "82hm", n_clones = 6, linesize = .4)

