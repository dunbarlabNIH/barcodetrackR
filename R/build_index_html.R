#' Build html
#'
#' Build html for vignette to index.html in docs
#' @param target the vignette to build
#' @param output the target for the vignette output
#'
#' @return Writes the vignette to docs/index.html. Only for internal use (get outta here!).
#'
# "
build_index_html <- function(target = "vignettes/Introduction_to_barcodetrackR.Rmd", output = "index.html") {
    rmarkdown::render(input = target, output_file = output, output_dir = "docs")
}
