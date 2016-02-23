#' Launch Multimonkey App
#'
#' Launches the Shiny Multimonkey App.
#'
#'@return Page launching the Shiny Barcode App
#'@examples
#'launchAppII()
#'@export
launchAppII <- function(){
  shiny::runApp(system.file("multimonkey_app", package = "barcodetrackR"))
}
