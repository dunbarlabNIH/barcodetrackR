% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_barcode_threshold.R
\name{estimate_barcode_threshold}
\alias{estimate_barcode_threshold}
\title{Estimate Barcode Threshold}
\usage{
estimate_barcode_threshold(
  capture_efficiency = NULL,
  population_size,
  proportion_labeled,
  confidence_level = 0.95,
  verbose = TRUE
)
}
\arguments{
\item{capture_efficiency}{Numeric. The capture efficiency of the clonal tracking method to detect a given clone. Must be between 0 and 1. See the description for details on how to estimate this value for a given experiment.}

\item{population_size}{Numeric. The total number of cells/genomes within each sample analyzed in the clonal tracking study. This is an experimentally determined value.}

\item{proportion_labeled}{Numeric. The proportion of the `population_size` which is genetically modified or contains a clonal tracking index. This is an experimentally determined value.}

\item{confidence_level}{Numeric. The confidence level for estimatig the minimum abundance threshold. Must be between 0 and 1. Default is 0.95 for 95 percent confidence that a clone with proportion `relative_threshold` will be detected. Increasing this parameter closer to one will result in a more stringent abundance threshold and decreasing this parameter will result in a more permissive abundance threshold.}

\item{verbose}{Logical. Whether to print the calculated threshold.}
}
\value{
Returns a single numeric `relative_threshold` describing the proportional abundance above which clones can be considered reliable given the provided capture efficiency and labeled population size. Pass this value into the function `threshold_SE` to threshold an existing SummarizedExperiment object or the function `create_SE` to threshold a SummarizedExperiment object upon creation from dataframes of counts and metadata.
}
\description{
Estimates an appropriate minimum abundance threshold for reliably detected barcodes in a clonal tracking dataset. \cr \cr
For a specified capture efficiency C, the minimum clone size N that we can expect to detect with confidence level P is calculated from: \cr `P = 1 - (1 - C)^(N)` \cr \cr
The proportional abundance of a clonal tag of size N is \cr `N / (T * F)` \cr where T is the total population size of cells or genomes and F is the frequency or proportion of the total population which is labeled or genetically modified with the clonal tag. \cr \cr
The population size and proportion labeled must be determined experimentally. The capture efficiency should be estimated for a given clonal tracking technique by simulating the barcode retrieval process in silico and finding the capture efficiency which leads to a total # of detected barcodes matching the experimentally determined number. Adair et al `(PMID: 32355868)` performed this analysis for viral integration site analysis and DNA barcode sequencing and determined good estimates for the capture efficiencies of these two technologies to be 0.05 and 0.4 respectively.
}
\examples{
estimate_barcode_threshold(
    capture_efficiency = 0.4,
    population_size = 500000,
    proportion_labeled = 0.3,
    confidence_level = 0.95,
    verbose = TRUE
)
}
