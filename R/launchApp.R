#' Launch Barcode App
#'
#' Launches the Shiny Barcode App.
#'
#'@return Page launching the Shiny Barcode App
#'@examples
#'launchApp()
#'@export
launchApp <- function(){
  shiny::runGitHub("d93espinoza/barcode_app")
}
