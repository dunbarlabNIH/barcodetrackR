#download the raw data from https://www.bbmt.org/article/S1083-8791(19)30566-X/abstract

require(dplyr)
require(magrittr)

raw_belderbos_data <- readxl::read_excel("inst/sample_data/BelderbosME_etal/mmc2.xlsx", sheet = "C21")
colnames(raw_belderbos_data) <- gsub("BM", "BM_", colnames(raw_belderbos_data))
raw_belderbos_data <- raw_belderbos_data[,!grepl("_sum", colnames(raw_belderbos_data))]
raw_belderbos_data <- as.data.frame(raw_belderbos_data)
rownames(raw_belderbos_data) <- raw_belderbos_data$ID
raw_belderbos_data$ID <- NULL
raw_belderbos_data <- raw_belderbos_data[,colSums(raw_belderbos_data) > 100]

belderbos_metadata <- data.frame(SAMPLENAME = colnames(raw_belderbos_data),
                                 weeks = c(9,14,20,22,22,22,rep("sac", 17)),
                                 celltype = c(rep("bulk", 4), "B", "T", "bulk",
                                              "B", "T", "B", "T", "G", "bulk",
                                              "B", "T", "G", "bulk", "B", "T", "G",
                                              "bulk", "B", "bulk"),
                                 organ = c(rep("PB", 6), rep("BM", 16), "Spleen"))

write.table(file = "inst/sample_data/BelderbosME_etal/count_matrix_mouse_C21.txt", raw_belderbos_data, sep = '\t', quote = FALSE)
write.table(file = "inst/sample_data/BelderbosME_etal/metadata_mouse_C21.txt", belderbos_metadata, sep = '\t', quote = FALSE)

