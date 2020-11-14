
shinyUI(
  fluidPage(
    titlePanel("Barcode App"),
    a(href="https://github.com/dunbarlabNIH/barcodetrackR", "Link to GitHub repository (source code, sample data, etc.)"),
    br(),
    a(href="http://dunbarlabNIH.github.io/barcodetrackR", "Link to vignette (learn how to use barcodetrackR)"),
    hr(),
    fluidRow(
      column(4,
             h4("1. Upload data"),
             div(style="display:inline-block; width:100%;",fileInput('file1', 'Tabular data with rows as observations, columns as samples'))
      ),
      column(4,
             h4("2. Upload metadata"),
             div(style="display:inline-block; width:100%;",fileInput('file2', 'Tabular metadata with SAMPLENAME column matching columns of data'))
      ),
      column(4,
             h4("3. Set threshold"),
             div(style="display:inline-block; width:100%;",numericInput("thresholdvalue", NULL, value = 0.0000)),
             div(style="display:inline-block; width:100%;",uiOutput("thresholdPanel")),
             div(style="display:inline-block; width:100%;", uiOutput("thresholdInfo"))
      ),
      #  div(style="display:inline-block; width:1%;", h5("")),
      #  div(style="display:inline-block; width:40%;",fileInput('file2', '2. Upload metadata (must have SAMPLENAME column matching columns of data)')),
      #  div(style="display:inline-block; width:1%;", h5("")),
      # # div(style="display:inline-block; width:15%;",fileInput('file3', '3. Upload readme (optional)')),
      # # div(style="display:inline-block; width:1%;", h5("")),
      #  div(style="display:inline-block; width:15%;",numericInput("thresholdvalue", "3. Set threshold", value = 0.0000)), #, div(style="display:inline-block; width:1%;", h5("")),
      #  div(style="display:inline-block; width:1%;", h5("")),
      #  div(style="display:inline-block; width:15%;",uiOutput("thresholdPanel")),
      #  div(style="display:inline-block; width:10%;", h5("")),
      #  div(style="display:inline-block; width:15%;", uiOutput("thresholdInfo")),
      # style = "padding: 10px;"
      
    ),
    uiOutput("Panel")
  )
)









