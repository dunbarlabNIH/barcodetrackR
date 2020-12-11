shinyUI(
  fluidPage(
    titlePanel("barcodetrackR Shiny app"),
    br(),
    tabsetPanel(id="Panel", type="tabs",
                tabPanel(title = "Upload Data",
                         fluidRow(column(4,
                                         wellPanel(
                                           h3("Option 1: Upload your data"),
                                           br(),
                                           br(),
                                           div(style="display:inline-block; width:100%;",fileInput('file1', 'Tabular data with rows as observations, columns as samples.')),
                                           div(style="display:inline-block; width:100%;",fileInput('file2', 'Tabular metadata with SAMPLENAME column matching columns of data.')),
                                           div(style="display:inline-block; width:100%;",
                                               numericInput("thresholdvalue",'Keep barcodes with proportion > threshold in at least one sample.',
                                                            value = 0.0000, min = 0, max = 1, step = 0.001)),
                                           div(style="display:inline-block; width:100%;",uiOutput("thresholdPanel")),
                                           div(style="display:inline-block; width:100%;", uiOutput("thresholdInfo")))),
                                  column(4,
                                         wellPanel(
                                           h3("Option 2: Load sample data"),
                                           br(),
                                           br(),
                                           actionButton("samplebutton", "Load Sample Data", width="100%"),
                                           br(),
                                           br(),
                                           br(),
                                           h5("What is the Sample Data?"),
                                           tags$em("Sample Data is a subset of barcodes and samples from the larger dataset ",
                                                   tags$a(href="https://github.com/dunbarlabNIH/barcodetrackR/tree/master/inst/sample_data/app_sample_data", "here")
                                                   )
                                           )
                                         ),

                                  column(4,
                                         wellPanel(
                                           h3("Links to help"),
                                           br(),
                                           tags$li(tags$a(href="http://dunbarlabNIH.github.io/barcodetrackR", "Link to vignette (learn how to use barcodetrackR).")),
                                           br(),
                                           tags$li(tags$a(href="https://github.com/dunbarlabNIH/barcodetrackR/tree/master/inst/sample_data/app_sample_data", "Link to format of example data.")),
                                           br(),
                                           tags$li(tags$a(href="https://github.com/dunbarlabNIH/barcodetrackR", "Link to GitHub repository (source code, sample data, issues, etc.).")),
                                           br(),
                                           br()))))
    )
  )

)
