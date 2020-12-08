
shinyUI(
  fluidPage(
    titlePanel("Barcode App"),
    a(href="https://github.com/dunbarlabNIH/barcodetrackR", "Link to GitHub repository (source code, sample data, issues, etc.)."),
    br(),
    a(href="http://dunbarlabNIH.github.io/barcodetrackR", "Link to vignette (learn how to use barcodetrackR)."),
    br(),
    hr(),
    fluidRow(
      column(4,
             h4("1. Upload data"),
             div(style="display:inline-block; width:100%;",fileInput('file1', 'Tabular data with rows as observations, columns as samples.')),
             div(style="display:inline-block; width:50%;",actionButton(inputId="loadSampleData", label="Load Sample Data")),
             uiOutput("SampleDataText")
      ),
      column(4,
             h4("2. Upload metadata"),
             div(style="display:inline-block; width:100%;",fileInput('file2', 'Tabular metadata with SAMPLENAME column matching columns of data.')),
             div(style="display:inline-block; width:50%;",actionButton(inputId="loadSampleMetadata", label="Load Sample Metadata")),
             uiOutput("SampleMetadataText")
      ),
      column(4,
             h4("3. Set threshold"),
             div(style="display:inline-block; width:100%;",numericInput("thresholdvalue",'Keep barcodes with proportion > threshold in at least one sample.', value = 0.0000, min = 0, max = 1, step = 0.001)),
             div(style="display:inline-block; width:100%;",uiOutput("thresholdPanel")),
             div(style="display:inline-block; width:100%;", uiOutput("thresholdInfo"))
      ),
    ),
    a(href="https://github.com/dunbarlabNIH/barcodetrackR/tree/master/inst/sample_data/app_sample_data", "Link to description of sample data."),
    hr(),
    uiOutput("Panel")
  )
)









