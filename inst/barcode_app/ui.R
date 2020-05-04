
shinyUI(
  fluidPage(
    titlePanel("Barcode App"),
    a(href="http://github.com/d93espinoza/barcodetrackR", "Link to GitHub repository (source code, sample data, etc.)"),
    h5(""),
    fluidRow(
      column(12,
             wellPanel(
               div(style="display:inline-block; width:40%;",fileInput('file1', '1. Upload data (rows as observations, columns as samples)')),
               div(style="display:inline-block; width:1%;", h5("")),
               div(style="display:inline-block; width:40%;",fileInput('file2', '2. Upload metadata (must have SAMPLENAME column matching columns of data)')),
               div(style="display:inline-block; width:1%;", h5("")),
              # div(style="display:inline-block; width:15%;",fileInput('file3', '3. Upload readme (optional)')),
              # div(style="display:inline-block; width:1%;", h5("")),
               div(style="display:inline-block; width:15%;",numericInput("thresholdvalue", "4. Set threshold", value = 0.0000)), div(style="display:inline-block; width:1%;", h5("")),
               div(style="display:inline-block; width:1%;", h5("")),
               div(style="display:inline-block; width:15%;",uiOutput("thresholdPanel")),
               div(style="display:inline-block; width:1%;", h5("")),
               div(style="display:inline-block; width:15%;", uiOutput("thresholdInfo")),
               style = "padding: 10px;"
             )
      )
    ),
    uiOutput("Panel")
  )
)









