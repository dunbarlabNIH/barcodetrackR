# barcodetrackR

This is an R package originally built for the use of cellular barcoding experiments performed at the Dunbar Lab at the National Institutes of Health. It includes a Shiny app for ease of use.

Code by Diego A. Espinoza and Samson J. Koelle.

## Installation

To install this package, try the following commands in R:
```
install.packages("devtools") #or skip, if you already have it
devtools::install_github("d93espinoza/barcodetrackR") #installs package from GitHub
```

## Running the app

After installation, use the following code to run the associated app:
```
barcodetrackR::launchApp()
```

Outfile must a tab-delimited .txt with barcodes as rows, samples as column headers and reads populating the sample matrix.
Below is the needed format:

|   | File 1 | File 2 | File 3 | File 4 | ... |
| ------------- | ------------- | ------------- | ------------- | ------------- | ----- |
| Barcode 1 | 1000 | 934 | 955 | 20 | ... |
| Barcode 2 | 450 | 90004 | 0 | 0 |... |
| Barcode 3  | 300 | 5001 | 95 | 100000 |... |
| ...  | ... | ... | ... | ... |... |



Keyfile must be a tab-delimited .txt file that must include the columns "FILENAME" and "GIVENNAME" as so (may include other columns):

| FILENAME | GIVENNAME |
| ------------- | ------------- |
| File 1 | Sample name 1 |
| File 2 | Sample name 2 | 
| File 3 | Sample name 3 |
| File 4  | Sample name 4 |
| ...  | ... |

The optional README must be a tab-delmited .txt. Lines may be commented out anywhere on the .txt and will not be read.
The README needs to have a FILENAME, MAPPED, and READS column.
Below is the example format:

|FILENAME | MAPPED | READS | ... |
| ------------- | ------------- | ------------- | ------------- |
| File 1 | 3000456 | 4000000 | ... |
| File 2 | 500000 | 3500000||... |
| File 3  | 100893 | 400000 | ... |
| ...  | ... | ... | ... |


## Included Functions

The following R functions are included in this package:

```
barcodecount()
#counts the number of barcodes across samples (unique, cumulative, new)

BBHM()
#shows the emergence of barcodes over time in a heatmap

barcode_ggheatmap()
#displays heatmap of the top N clones in selected samples

clonaldiversity()
#calcuate diversity indices for samples across time

cor_plot()
#wrapper for corrplot from corrplot package, with normalizations for barcode data

diversity_plotter()
#plot the diversity of samples across time

gettopindices()
#extracts the top N indices for a data frame per column, eliminating repeats

launchApp()
#launch the shiny app that eases the use of the majority these functions

merger()
#merges list of matrices by rowname, including all entries

radartopclones()
#radar plot of top N clones in a sample across samples

richness_plotter()
#plots the barcode counts of samples across time

threshold()
#eliminates barcodes that aren't present X times in at least one sample

treemap()
#tree map of a given sample
```
