# diego espinoza
# 12/16/2019

# Purpose: to load barcode data and associated metadata into SummarizedExperiment
# to call our different visualization functions on data stored in this object

require(magrittr)

# Load barcode and meta data
zj31_heatmap_samples <- c("zj316mPBT28cyclesCD75R.fastq", "zj31_9hm_PB_T_CD3p_CD56n_R23_COMBINEDCD99CD100.fastq",
                          "zj31_12m_T_CD_103_R3_TTAGGC_L004_R1_001.fastq", "zj31_20m_T_CD3p_20160328_pp_50trim_CD124_R2.fastq",
                          "zj316mPBB28cyclesCD75R.fastq", "zj31_9hm_PB_B_R24_COMBINEDCD99CD100.fastq",
                          "zj31_12m_B_CD_103_R4_TGACCA_L004_R1_001.fastq", "zj31_20m_B_20160328_pp_50trim_CD124_R3.fastq",
                          "zj316mPBGr28cyclesCD75R.fastq", "zj31_9hm_PB_Gr_CD33p_GFPnp_R22_COMBINEDCD99CD100.fastq",
                          "zj31_12m_Gr_CD33p_CD_103_R2_CGATGT_L004_R1_001.fastq", "zj31_20m_Gr_20160328_pp_50trim_CD124_R1.fastq",
                          "zj316mCD56spNK28cyclesCD76R.fastq", "zj31_9hm_PB_CD56p_CD100_R27_ATTCCT_L007_R1_001.fastq",
                          "zj31_12m_20150723_frozen_CD56sp_2050817SORT_CD_105_R4_TGACCA_L007_R1_001.fastq",
                          "zj31_20m_NK_CD56sp_20160328_pp_50trim_CD124_R6.fastq", "zj316mCD16sp28cyclesCD76R.fastq",
                          "zj31_9hm_PB_CD16p_R26_COMBINEDCD99CD100.fastq", "zj31_12m_CD16sp_CD_103_R6_GCCAAT_L004_R1_001.fastq",
                          "zj31_20m_NK_CD16sp_20160328_pp_50trim_CD124_R5.fastq")
barcode.file <- read.delim("~/Documents/dunbar_lab/projects/joy_NK/manuscript/SI_review/final_outfiles/ZJ31_combined_20180828_NKpaper.txt", row.names = 1) %>% .[,zj31_heatmap_samples]
zj31_fig2_heatmap_names <- paste(c("6m", "9.5m", "12m", "20m"), rep(c("T", "B", "Gr", "NK CD56+ CD16-", "NK CD56- CD16+"), each = 4))
zj31_timepoints <- rep(c(6, 9.5, 12, 20), times = 5)
zj31_celltypes <- rep(c("T", "B", "Gr", "NK 56", "NK 16"), each = 4)

meta.file <- data.frame(SAMPLENAME = zj31_heatmap_samples, mynames = zj31_fig2_heatmap_names, Lineage = zj31_celltypes, Month = zj31_timepoints)

# Create summarized experiment
SE <- create_SE(your_data = barcode.file, meta_data = meta.file)

# Heat map of just the B cells
barcodetrackR::barcode_ggheatmap_2(your_data = SummarizedExperiment::assay(se)[,se$Cell_type == "B"],
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

