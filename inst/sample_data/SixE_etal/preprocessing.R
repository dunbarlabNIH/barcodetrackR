#download the raw data from https://github.com/BushmanLab/HSC_diversity/blob/master/data/intSites.mergedSamples.collapsed.csv.gz which is the data repository for https://ashpublications.org/blood/article-abstract/doi/10.1182/blood.2019002350/441042/Clonal-tracking-in-gene-therapy-patients-reveals-a?redirectedFrom=fulltext


raw_six_data <- read.csv("inst/sample_data/SixE_etal/intSites.mergedSamples.collapsed.csv.gz")

raw_six_data %>%
  dplyr::filter(patient == "WAS5") %>%
  dplyr::mutate(new_col_id = paste0("m", timePointMonths, "_", cellType)) %>%
  tidyr::pivot_wider(id_cols = posid, values_from = reads, names_from = new_col_id) %>%
  as.data.frame -> raw_six_data_reads

rownames(raw_six_data_reads) <- raw_six_data_reads$posid
raw_six_data_reads$posid <- NULL
raw_six_data_reads[is.na(raw_six_data_reads)] <- 0

raw_six_data %>%
  dplyr::filter(patient == "WAS5") %>%
  dplyr::mutate(new_col_id = paste0("m", timePointMonths, "_", cellType)) %>%
  tidyr::pivot_wider(id_cols = posid, values_from = estAbund, names_from = new_col_id) %>%
  as.data.frame -> raw_six_data_estabundance

rownames(raw_six_data_estabundance) <- raw_six_data_estabundance$posid
raw_six_data_estabundance$posid <- NULL
raw_six_data_estabundance[is.na(raw_six_data_estabundance)] <- 0


colname_order <- c("m13_TCELLS", "m13_BCELLS", "m13_NKCELLS", "m13_GRANULOCYTES", "m13_MONOCYTES",
                   "m36_TCELLS", "m36_BCELLS",  "m36_NKCELLS",  "m36_GRANULOCYTES", "m36_MONOCYTES",
                   "m43_TCELLS", "m43_BCELLS", "m43_NKCELLS", "m43_GRANULOCYTES", "m43_MONOCYTES",
                   "m55_TCELLS", "m55_BCELLS",  "m55_NKCELLS", "m55_GRANULOCYTES", "m55_MONOCYTES")
raw_six_data_estabundance <- raw_six_data_estabundance[,colname_order]
raw_six_data_reads <- raw_six_data_reads[,colname_order]


WAS5_metadata <- data.frame(SAMPLENAME = colnames(raw_six_data_estabundance),
                            months = rep(c(13, 36, 43, 55), each = 5),
                            celltype = rep(c("T", "B", "NK", "Gr", "Mo"), times = 4))


all(rownames(raw_six_data_estabundance) == rownames(raw_six_data_reads))

write.table(file = "inst/sample_data/SixE_etal/WAS5_reads.txt", raw_six_data_reads, sep = '\t', quote = FALSE)
write.table(file = "inst/sample_data/SixE_etal/WAS5_estabundance.txt", raw_six_data_estabundance, sep = '\t', quote = FALSE)
write.table(file = "inst/sample_data/SixE_etal/WAS5_metadata.txt", WAS5_metadata, sep = '\t', quote =  FALSE, row.names = FALSE)
