# Ryland Mortlock
# 12/13/2019

# Purpose: to load barcode data and associated metadata into SummarizedExperiment
          # to call our different visualization functions on data stored in this object

# Load and threshold data
file <- '/Users/mortlockrd/Desktop/GitHub/barcodetrackR/data/ZG66_simple_data.txt'
meta_file <- '/Users/mortlockrd/Desktop/GitHub/barcodetrackR/data/ZG66_simple_metadata.txt'
bc.df <- barcodetrackR::threshold(read.delim(file, row.names = 1),thresh = 0.0005)

# Load metadata
meta.df <- read.delim(meta_file)

# Create summarized experiment
bc.mat <- as.matrix(bc.df)
se <- SummarizedExperiment::SummarizedExperiment(assays = list(bc.mat = bc.mat),colData=meta.df)

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

# Clone bar plot
source('clone_bar_plot.R')
clone_bar_plot(your_SE = se, filter_by = "Cell_type", filter_selection = "T", plot_by = "Timepoint", n_clones = 20, linesize = .4)
