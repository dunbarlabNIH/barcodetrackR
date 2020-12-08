#download the raw data from https://doi.org/10.5281/zenodo.1256169 which is the data repository for https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0580-z

## T Cell Receptor (TCR) Sequencing Data
load('your_file_path/clarke_tcr_data.rds') # Downloaded from zenodo link above
# The TCR data is too large to include in the Github repository, it can be downloaded from the link above and saved locally. 

# Counts data
tcr[["seqs"]] %>%
  mutate(identifier = paste(accn, replicate, sep = "_")) %>%
  tidyr::pivot_wider(id_cols = nucleotide, values_from = count, names_from = identifier) %>%
  as.data.frame -> clarke_TCR_reads

rownames(clarke_TCR_reads) <- clarke_TCR_reads$nucleotide
clarke_TCR_reads$nucleotide <- NULL
clarke_TCR_reads[is.na(clarke_TCR_reads)] <- 0

# Metadata keeping only patient 7
tcr[["mdata"]] %>%
  dplyr::filter(patient == "SCID00007") %>%
  mutate(SAMPLENAME = paste(accn, replicate, sep = "_")) %>%
  arrange(nmonths, replicate) %>%
  select(SAMPLENAME, everything()) -> clarke_TCR_metadata_p7

# Keep just patient 7 data in read counts
clarke_TCR_reads_p7 <- clarke_TCR_reads[,colnames(clarke_TCR_reads) %in% clarke_TCR_metadata_p7$SAMPLENAME]
clarke_TCR_reads_p7 <- clarke_TCR_reads_p7[rowSums(clarke_TCR_reads_p7)>0,]

# Match order from metadata
clarke_TCR_reads_p7 <- clarke_TCR_reads_p7[,match(colnames(clarke_TCR_reads_p7),clarke_TCR_metadata_p7$SAMPLENAME)]

# Write to files
write.table(file = "inst/sample_data/ClarkeEL_etal/tcr_reads_p00007.txt", clarke_TCR_reads_p7, sep = '\t', quote = FALSE)
write.table(file = "inst/sample_data/ClarkeEL_etal/tcr_metadata_p00007.txt", clarke_TCR_metadata_p7, sep = '\t', quote =  FALSE)


## Viral Integration Site (VIS) Data
raw_clarke_data <- read.csv("inst/sample_data/ClarkeEL_etal/intSiteData.csv")

raw_clarke_data %>%
  dplyr::filter(patient == "p00007") %>%
  mutate(identifier = paste(chromosome, position, strand, sep = ".")) %>%
  tidyr::pivot_wider(id_cols = identifier, values_from = reads, names_from = timePointMonths) %>%
  as.data.frame -> raw_clarke_data_reads

rownames(raw_clarke_data_reads) <- raw_clarke_data_reads$identifier
raw_clarke_data_reads$identifier <- NULL
colnames(raw_clarke_data_reads) <- paste0(colnames(raw_clarke_data_reads), "_months")
raw_clarke_data_reads[is.na(raw_clarke_data_reads)] <- 0

raw_clarke_data %>%
  dplyr::filter(patient == "p00007") %>%
  mutate(identifier = paste(chromosome, position, strand, sep = ".")) %>%
  tidyr::pivot_wider(id_cols = identifier, values_from = inferredCells, names_from = timePointMonths) %>%
  as.data.frame -> raw_clarke_data_cells

rownames(raw_clarke_data_cells) <- raw_clarke_data_cells$identifier
raw_clarke_data_cells$identifier <- NULL
colnames(raw_clarke_data_cells) <- paste0(colnames(raw_clarke_data_cells), "_months")
raw_clarke_data_cells[is.na(raw_clarke_data_cells)] <- 0

p00007.metadata <- data.frame(SAMPLENAME = c("6_months", "12_months", "18_months", "24_months"),
                              months = c(6, 12, 18, 24))

all(rownames(raw_clarke_data_cells) == rownames(raw_clarke_data_reads))

write.table(file = "inst/sample_data/ClarkeEL_etal/p00007_reads.txt", raw_clarke_data_reads, sep = '\t', quote = FALSE)
write.table(file = "inst/sample_data/ClarkeEL_etal/p00007_cells.txt", raw_clarke_data_cells, sep = '\t', quote = FALSE)
write.table(file = "inst/sample_data/ClarkeEL_etal/p00007_metadata.txt", p00007.metadata, sep = '\t', quote =  FALSE)
