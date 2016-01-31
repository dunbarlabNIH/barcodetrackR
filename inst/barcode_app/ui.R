
shinyUI(
  fluidPage(
    titlePanel("Barcode App"),
    a(href="http://github.com/d93espinoza/", "Link to GitHub repository (source code, sample data, etc.)"),
    h5(""),
    br(),
    
    fluidRow(
      column(3,
             wellPanel(
               width = 12,
               fileInput('file1', '1. Choose Barcode OutFile')
             ) 
      ),
      column(3,
             wellPanel(
               fileInput('file2', '2. Choose Your KeyFile')
               
             )
      ), 
      
      column(3,
             wellPanel(
               numericInput("thresholdvalue", "3. Threshold Value to Apply", value = 2000)
             )
      ),
      
      uiOutput("thresholdPanel")
      
    ),#end of fluidRow
    
    
    fluidRow(
      
      uiOutput("thresholdInfo")
    
    ),
    
    uiOutput("Panel")
  )
)









