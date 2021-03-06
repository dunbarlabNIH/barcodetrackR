% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rank_abundance_stat_test.R
\name{rank_abundance_stat_test}
\alias{rank_abundance_stat_test}
\title{Rank Abundance Statistical Test}
\usage{
rank_abundance_stat_test(your_SE, statistical_test = "ks")
}
\arguments{
\item{your_SE}{Summarized Experiment object containing clonal tracking data as created by the barcodetrackR `create_SE` function.}

\item{statistical_test}{The statistical test used to compare distributions. For now, the only implemented test is the Kolmogorov-Smirnov test.}
}
\value{
Returns a list containing two dataframes \cr [["D_statistic"]] is a dataframe containing pairwise D-statistics between each pair of samples in your_SE. The D statistic represents the maximal difference between the two rank abundance distributions. \cr [["p_value]] A dataframe containing the p-value computed by the KS test for each pair of samples. The null hypothesis is that the two rank-abundance profiles come from the same distribution.
}
\description{
Carries out a specific instance of statistical testing relevant to clonal tracking experiments. For the provided SummarizedExperiment,
compare the rank-abundance distribution which is described by the increase in cumulative abundance within that sample as barcode abundances are added, starting with the most abundant barcode.
The two-sided Kolmogorov-Smirnov statistical test is carried out comparing each pair of samples using the R function ks.test:https://www.rdocumentation.org/packages/dgof/versions/1.2/topics/ks.test
Note that this test compares rank-abundance distribution regardless of whether the samples share the same barcodes or lineage tracing elements. The test could be employed on two samples with no barcode sequence overlap, simply to compare whether the rank abundance distribution of barcodes is drawn from the same distribution.
}
\examples{
data(wu_subset)
rank_abundance_stat_test(your_SE = wu_subset, statistical_test = "ks")
}
