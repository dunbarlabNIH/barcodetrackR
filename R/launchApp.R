#' Launch Barcode App
#'
#' Launches the Shiny Barcode App.
#'
#' @param x NULL
#' @return Page launching the Shiny Barcode App
#' @examples
#' \dontrun{
#' launchApp()
#' }
#' @export
launchApp <- function(x = NULL) {
    shiny::runApp(system.file("barcode_app", package = "barcodetrackR"))
}
