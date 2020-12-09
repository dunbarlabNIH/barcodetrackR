#' Build html
#'
#' Build html for vignette to index.html in docs
#'
#"
build_index_html <- function(target = "vignettes/Introduction_to_barcodetrackR.Rmd", output = "index.html"){
  rmarkdown::render(input = target, output_file = output, output_dir = "docs")
}
