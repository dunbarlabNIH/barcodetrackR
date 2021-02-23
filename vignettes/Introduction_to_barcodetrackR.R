## ----setup, echo = FALSE------------------------------------------------------
#knitr::opts_chunk$set(width = 100)

## ----installation, eval = FALSE-----------------------------------------------
#  devtools::install_github("dunbarlabNIH/barcodetrackR")

## ----load packages required, eval = TRUE, warning=FALSE, message=FALSE--------

require("magrittr")
require("barcodetrackR")
require("SummarizedExperiment")


## ----load data, eval = TRUE---------------------------------------------------

system.file("sample_data/WuC_etal/monkey_ZJ31.txt", package = "barcodetrackR") %>%
  read.delim(row.names = 1) -> wu_dataframe

system.file("sample_data/WuC_etal/monkey_ZJ31_metadata.txt", package = "barcodetrackR") %>%
  read.delim() -> wu_metadata

wu_SE <- create_SE(your_data = wu_dataframe,
                   meta_data = wu_metadata,
                   threshold = 0)

system.file("sample_data/KoelleSJ_etal/ZH33_reads.txt", package = "barcodetrackR") %>%
  read.delim(row.names = 1) -> koelle_dataframe

system.file("sample_data/KoelleSJ_etal/ZH33_metadata.txt", package = "barcodetrackR") %>%
  read.delim() -> koelle_metadata

koelle_SE <- create_SE(your_data = koelle_dataframe,
                       meta_data = koelle_metadata,
                       threshold = 0)

system.file("sample_data/BelderbosME_etal/count_matrix_mouse_C21.txt", package = "barcodetrackR") %>%
  read.delim(row.names = 1) -> belderbos_dataframe

system.file("sample_data/BelderbosME_etal/metadata_mouse_C21.txt", package = "barcodetrackR") %>%
  read.delim() -> belderbos_metadata

belderbos_metadata$weeks <- factor(belderbos_metadata$weeks, levels = c("9", "14", "20", "22", "sac"))

belderbos_SE <- create_SE(your_data = belderbos_dataframe,
                          meta_data = belderbos_metadata,
                          threshold = 0)



system.file("sample_data/SixE_etal/WAS5_reads.txt", package = "barcodetrackR") %>%
  read.delim(row.names = 1) -> six_dataframe

system.file("sample_data/SixE_etal/WAS5_metadata.txt", package = "barcodetrackR") %>%
  read.delim() -> six_metadata

six_SE <- create_SE(your_data = six_dataframe,
                    meta_data = six_metadata,
                    threshold = 0)


## ----list assays, eval = TRUE-------------------------------------------------
assays(six_SE)

## ----estimate threshold, eval = TRUE------------------------------------------
my_thresh <- estimate_barcode_threshold(capture_efficiency = 0.4, population_size = 500000, proportion_labeled = 0.3,
                                        confidence_level = 0.95, verbose = TRUE)

## ----apply thresholding, eval = TRUE------------------------------------------
wu_thresh <- threshold_SE(your_SE = wu_SE, threshold_value = 0.005, threshold_type = "relative",
                          verbose = TRUE)

## ----launch_app, eval=FALSE---------------------------------------------------
#  barcodetrackR::launchApp()

## ----scatter_plot, eval = TRUE, fig.width = 7, fig.height=4-------------------
Gr_B_20 <- c("ZJ31_20m_Gr", "ZJ31_20m_B")
Gr_T_20 <- c("ZJ31_20m_Gr", "ZJ31_20m_T")
wu_scatterplot_1 <- scatter_plot(wu_SE[,Gr_B_20], your_title = "Gr vs B", assay = "proportions")
wu_scatterplot_2 <- scatter_plot(wu_SE[,Gr_T_20], your_title = "Gr vs T", assay = "proportions")
cowplot::plot_grid(wu_scatterplot_1, wu_scatterplot_2, ncol = 2)

## ----cor_plot 1, eval = TRUE, fig.width = 6, fig.height = 5-------------------
wu_cor_plot_sample_selection <- colData(wu_SE)$SAMPLENAME[1:20]
cor_plot(wu_SE[,wu_cor_plot_sample_selection],
         method_corr = "pearson",
         plot_type = "color",
         assay = "proportions")

## ----cor_plot 2, eval = TRUE, fig.width = 6, fig.height = 5-------------------
wu_cor_plot_sample_selection <- colData(wu_SE)$SAMPLENAME[1:20]
cor_plot(wu_SE[,wu_cor_plot_sample_selection],
         method_corr = "pearson",
         plot_type = "color",
         assay = "logs")

## ----cor_plot 3, eval = TRUE--------------------------------------------------
cor_plot(wu_SE[,wu_cor_plot_sample_selection],
         method_corr = "pearson",
         assay = "proportions",
         return_table = TRUE) %>% head

## ----cor_plot 4, eval = TRUE, fig.width = 4, fig.height = 3, warning=FALSE, message=FALSE----
belderbos_wk_samples_PB <- c("wk9_U", "wk14_U", "wk20_U", "wk22_U")
cor_plot(your_SE = belderbos_SE[,belderbos_wk_samples_PB],
         method_corr = "pearson",
         plot_type = "number",
         assay = "proportions",
         number_size = 3)

## ----summary proxy, eval = TRUE-----------------------------------------------
summary(proxy::pr_DB)

## ----dist plot 1, eval = TRUE, fig.width = 7, fig.height = 5------------------
dist_plot(wu_SE[,wu_cor_plot_sample_selection],
          dist_method =  "manhattan",
          plot_type = "color",
          assay = "logs")

## ----dist plot 2, eval = TRUE, fig.width = 7, fig.height = 5------------------
dist_plot(wu_SE[,wu_cor_plot_sample_selection],
          dist_method =  "cosine",
          plot_type = "color",
          assay = "logs",
          color_pal = "Greens",
          cluster_tree = TRUE)

## ----barcode_ggheatmap wu_SE 1, echo=TRUE, warning = FALSE, fig.width = 12, fig.height=10----
barcode_ggheatmap(your_SE = wu_SE[,wu_cor_plot_sample_selection],
                  n_clones = 10,
                  label_size = 14,
                  cellnote_size = 4)

## ----barcode_ggheatmap six_SE, echo=TRUE, warning = FALSE, fig.width = 12, fig.height=10----
six_celltype_order <- c("m13_TCELLS", "m36_TCELLS", "m43_TCELLS", "m55_TCELLS",
                        "m13_BCELLS", "m36_BCELLS", "m43_BCELLS", "m55_BCELLS",
                        "m13_NKCELLS", "m36_NKCELLS", "m43_NKCELLS","m55_NKCELLS",
                        "m13_GRANULOCYTES", "m36_GRANULOCYTES", "m43_GRANULOCYTES",
                        "m55_GRANULOCYTES",
                        "m13_MONOCYTES", "m36_MONOCYTES", "m43_MONOCYTES","m55_MONOCYTES")

barcode_ggheatmap(your_SE = six_SE[,six_celltype_order],
                  n_clones = 5,
                  cellnote_assay = "stars",
                  cellnote_size = 6,
                  label_size = 14,
                  dendro = TRUE,
                  clusters = 4,
                  distance_method = "euclidean")

## ----barcode_ggheatmap belderbos_SE, echo=TRUE, warning = FALSE, fig.width = 7, fig.height = 5----
belderbos_wk_samples_PB <- c("wk9_U", "wk14_U", "wk20_U", "wk22_U")
barcode_ggheatmap(your_SE = belderbos_SE[,belderbos_wk_samples_PB],
                  n_clones = 10,
                  label_size = 18,
                  dendro = FALSE,
                  cellnote_size = 4,
                  cellnote_assay = "proportions")

## ----barcode_ggheatmap_stat wu_SE, echo=TRUE, warning = FALSE, fig.width = 10, fig.height=7----
wu_CD16_NK_order <- c("ZJ31_3m_NK_CD56n_CD16p", "ZJ31_4m_NK_CD56n_CD16p",
                      "ZJ31_6m_NK_CD56n_CD16p", "ZJ31_8.5m_NK_CD56n_CD16p",
                      "ZJ31_9.5m_NK_CD56n_CD16p", "ZJ31_12m_NK_CD56n_CD16p",
                      "ZJ31_17.5m_NK_CD56n_CD16p", "ZJ31_20m_NK_CD56n_CD16p")

barcode_ggheatmap_stat(your_SE = wu_SE[,wu_CD16_NK_order],
                       sample_size = rep(40000,length(wu_CD16_NK_order)),
                       stat_test = "chi-squared",
                       stat_option = "subsequent",
                       p_threshold = 0.05,
                       n_clones = 10,
                       cellnote_assay = "stars",
                       cellnote_size = 6)

## ----barcode_stat_test wu_SE, eval = TRUE, echo=TRUE, warning = FALSE---------

wu_CD16_NK_statistics <- barcode_stat_test(your_SE = wu_SE[,wu_CD16_NK_order],
                                           sample_size = rep(40000,length(wu_CD16_NK_order)),
                                           stat_test = "chi-squared",
                                           stat_option = "subsequent",
                                           bc_threshold = 0.01)

head(wu_CD16_NK_statistics$p_val[,4:5])

## ----binary_heat_map, echo=TRUE, warning = FALSE, eval = TRUE-----------------
barcode_binary_heatmap(your_SE = belderbos_SE[,belderbos_wk_samples_PB],
                       label_size = 12,
                       threshold = 0.01)

## ----clonal_contribution wu_SE, echo=TRUE, warning = FALSE--------------------
clonal_contribution(your_SE = wu_SE,
                    SAMPLENAME_choice = "ZJ31_20m_NK_CD56n_CD16p",
                    n_clones = 10,
                    graph_type = "line",
                    plot_over = "months",
                    filter_by = "celltype",
                    filter_selection = "NK_16",
                    plot_non_selected = FALSE)

