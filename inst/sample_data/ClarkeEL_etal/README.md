For using the Clarke et al sample data in the create_SE function or the barcodetrackR app, one can analyze:

* viral integration site counts using the count data file "p00007_reads.txt" and metadata file "p00007_metadata.txt"
* T Cell Receptor clonal counts using the count data file "tcr_reads_p00007.txt" and metadata file "tcr_metadata_p00007.txt""

The file "p00007_cells.txt" matches the format of the viral integration site counts but contains inferred number of cells.

The publicly available raw data from the published manuscript is available at https://zenodo.org/record/1256169#.X53beUJKiL4 <br /> The script preprocessing.R reads in this raw data and outputs counts and metadata files in the format amenable to barcodetrackR. Note that the raw data "intSiteData.csv"" from the prior mentioned link is included in the sample data repository but the TCR raw dataset mustbe downloaded locally because it is too large. 
