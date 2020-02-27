#' Launch Barcode App
#'
#' Launches the Shiny Barcode App.
#'
#'@return Page launching the Shiny Barcode App
#'@examples
#'launchApp()
#'@export
launchApp <- function(){
  shiny::runApp(system.file("barcode_app", package = "barcodetrackR"))
}
