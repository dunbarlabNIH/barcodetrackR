The sample data for the barcodetrackR app contains genetic barcoding data from the following published paper:
[Wu, Chuanfeng, et al. "Clonal expansion and compartmentalized maintenance of rhesus macaque NK cell subsets." Science Immunology (2018)](http://dx.doi.org/10.1126/sciimmunol.aat9781)

The file "sample_data_ZJ31.txt" contains the counts for approximately 35,000 genetic barcodes across the first 20 samples from a single animal which was autologously transplanted with lentivirally barcoded hematopoietic stem and progenitor cells in the Wu et al study. The counts represent barcode abundances quantified by high-throughput sequencing from DNA extracted from sorted lineage cells from the animal's peripheral blood. 

The file "sample_metadata_ZJ31.txt" describes the timepoint and cell type of the 20 samples. The cell types included are T (T cells), B (B cells), Gr (Granulocytes), NK_56 (CD56 NK cells), and NK_16 (CD16 NK cells). The timepoints included are 6, 9.5, 12, and 20 months post-transplant.

These files can be opened with a text editor or Excel to understand the structure of the input data for the barcodetrackR application. Note that the metadata file must contain a column named SAMPLENAME which matches with the column names of the counts data. 