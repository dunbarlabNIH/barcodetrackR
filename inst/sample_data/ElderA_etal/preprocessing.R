#download the raw data from GSEq49170
require(magrittr)
require(dplyr)
file_list <- list.files("GSE149170_RAW/")
count_list <- lapply(1:length(file_list), function(i){
  file_path <- paste0("GSE149170_RAW/", file_list[i])
  message("reading file ", file_list[i])
  new_name <- gsub(pattern = ".counts.xlsx", replacement = "", x= file_list[i])
  return(readxl::read_excel(file_path, col_names = c("barcode", new_name), range = readxl::cell_cols("A:B")))
})

full_matrix <- Reduce(function(...) dplyr::left_join(..., by='barcode'),count_list)
full_matrix[is.na(full_matrix)] <- 0


adf <- Biobase::phenoData(GEOquery::getGEO(filename = "GSE149170_series_matrix.txt", GSEMatrix = TRUE, getGPL = FALSE))
my_metadata <- adf@data

colnames(full_matrix) <- lapply(strsplit(colnames(full_matrix), split = "_"), function(x){x[1]}) %>% unlist
full_matrix <- tibble::column_to_rownames(full_matrix, var = "barcode")


mouse <- sapply(strsplit(my_metadata$title,"(", fixed = TRUE), `[`, 1) %>% gsub(pattern = " ", replacement = "", fixed = TRUE)
organ <- sapply(strsplit(my_metadata$title,"(", fixed = TRUE), `[`, 2) %>% gsub(pattern = ")", replacement = "", fixed = TRUE)
sample_type <- my_metadata$`sample type:ch1`
sample_origin <- my_metadata$`sample:ch1`
cell_dose <- my_metadata$`cells transplanted:ch1`
filtered_metadata <- data.frame(my_metadata$geo_accession, mouse, organ, sample_type, sample_origin, cell_dose)


filtered_metadata <- select(my_metadata, geo_accession, characteristics_ch1, characteristics_ch1.1, characteristics_ch1.2, characteristics_ch1.3, "sample type:ch1", "sample:ch1", "transduction:ch1")

write.table(file = "GSE148170_matrix.txt", full_matrix, sep = '\t', quote = FALSE, row.names = TRUE)
write.table(file = "GSE148170_metadata.txt", filtered_metadata, sep = '\t', quote = FALSE, row.names = FALSE)

my_metadata <- read.delim("GSE148170_metadata.txt")
my_counts <- read.delim("GSE148170_matrix.txt")

L4951_metadata <- my_metadata %>% filter(sample_origin == "L4951") %>% mutate(SAMPLENAME = my_metadata.geo_accession)
L4951_counts <- my_counts[,L4951_metadata$my_metadata.geo_accession]


short_organ <- L4951_metadata$organ %>% gsub("left", "L", .) %>% gsub("right", "R", .) %>% gsub(" part", "", .)
short_sample_type <- L4951_metadata$sample_type %>% gsub("primary transplant", "pri", .) %>% gsub("secondary of AE12 spleen", "sec AE12", .) %>% gsub("secondary of AE6 spleen", "sec AE6", .) %>% gsub("tertiary of AE32 spleen", "ter AE32", .) %>% gsub("tertiary of AE36 spleen", "ter AE36", .)

L4951_metadata$short_name <- paste0(L4951_metadata$mouse, " ", short_organ, " ", short_sample_type)

write.table(file = "L4951_count_matrix.txt", L4951_counts, sep = '\t', quote = FALSE, row.names = TRUE)
write.table(file = "L4951_metadata.txt", L4951_metadata, sep = '\t', quote = FALSE, row.names = FALSE)
