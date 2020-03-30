#download the raw data from https://github.com/dunbarlabNIH/R-code-and-tabular-data/tree/master/datacoderesults/data which is the data repository for https://ashpublications.org/blood/article-lookup/doi/10.1182/blood-2016-07-728691

koelle_metadata <- read.delim("inst/sample_data/KoelleSJ_etal/zh33keyfileforpub.txt", stringsAsFactors = FALSE) %>%
  dplyr::filter(!is.na(FILENAME) & CELLTYPE %in% c("T", "Gr", "Mono")) %>%
  dplyr::mutate(SAMPLENAME = paste0("ZH33_", gsub(" ", "_", GIVENNAME))) %>%
  dplyr::mutate(SAMPLENAME = gsub("_112514", "", SAMPLENAME)) %>%
  dplyr::mutate(MONTH = gsub("150806", "38", MONTH))

raw_koelle_data <- read.delim("inst/sample_data/KoelleSJ_etal/zh33_independent_100_forpub.txt")[,koelle_metadata$FILENAME]
raw_koelle_data <- raw_koelle_data[rowSums(raw_koelle_data) > 0,]
colnames(raw_koelle_data) <- plyr::mapvalues(colnames(raw_koelle_data),
                                             from = koelle_metadata$FILENAME,
                                             to = koelle_metadata$SAMPLENAME)

koelle_metadata %>%
  dplyr::rename(month = MONTH, celltype = CELLTYPE) %>%
  dplyr::select(SAMPLENAME, month, celltype) -> koelle_metadata


write.table(file = "inst/sample_data/KoelleSJ_etal/ZH33_reads.txt", raw_koelle_data, sep = '\t', quote = FALSE)
write.table(file = "inst/sample_data/KoelleSJ_etal/ZH33_metadata.txt", koelle_metadata, sep = '\t', quote = FALSE)
