#barcodetrackR

This is an R package originally built for the use of cellular barcoding experiments performed at the Dunbar Lab at the National Institutes of Health. It includes a Shiny app for ease of use.

Code by Diego A. Espinoza and Samson J. Koelle.

To install this package, try the following commands in R:
```
install.packages("devtools") #or skip, if you already have it
devtools::install_github("d93espinoza/barcodetrackR") #installs package from GitHub
```

After installation, use the following code to run the associated app:
```
barcodetrackR::launchApp()
```

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
