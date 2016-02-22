library(shinyBS)
shinyUI(
  fluidPage(
    titlePanel("Multimonkey App"),
    a(href="http://github.com/d93espinoza/", "Link to GitHub repository (source code, sample data, etc.)"),
    h5(""),

    #========================================================================================================================

    fluidRow(
      column(6,
             wellPanel(style="text-align:center",
               strong("1. Press button to upload monkey outfiles, keyfiles, and thresholds"),
               br(),
               actionButton("upload", "Upload Data")
             )
      ),
      column(6,
             wellPanel(style="text-align:center",
               strong("2. Press button to apply thresholds and begin analysis"),
               br(),
               actionButton("start_button", "Begin Analysis")

             ))
    ),


    shinyBS::bsModal("modalExample", "Upload Data", "upload", size = "large",
                     fluidRow(column(4,numericInput("nummonkeys", label = "Number of Monkeys", value = 1, min = 1))),
                     uiOutput("rowsofmonks")
    ),

    uiOutput("panel")


  )
)









