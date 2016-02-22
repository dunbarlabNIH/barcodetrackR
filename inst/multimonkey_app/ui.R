library(shinyBS)
shinyUI(
  fluidPage(
    titlePanel("Multimonkey App"),
    a(href="http://github.com/d93espinoza/", "Link to GitHub repository (source code, sample data, etc.)"),
    h5(""),

    fluidRow(
      column(3,
             wellPanel(
               strong("1. Press button to Upload Data"),
               actionButton("upload", "Upload Data")
             )
      ),
      column(3,
             wellPanel(style="height:97px",
                       numericInput("thresholdvalue", "2. Threshold Value to Apply", value = 2000)
             )
      ),
      column(3,
             wellPanel(
               strong("3. Press button to apply threshold"),
               br(),
               actionButton("threshold_button", "Apply Threshold")

             ))
    ),

    fluidRow(
      column(12,
             wellPanel(
               div(style="display:inline-block",
                   textInput("DS", "Monkey 1 Name")),
               div(style="display:inline-block",
                   fileInput("DSF", "DSFD")),
               div(style="display:inline-block;width=100%",
                   h5("Accuracy table:")))
      )

    ),

    shinyBS::bsModal("modalExample", "Upload Data", "upload", size = "large",
                     fluidRow(column(4,numericInput("nummonkeys", label = "Number of Monkeys", value = 1, min = 1))),
                     uiOutput("rowsofmonks")




    )

  )
)









