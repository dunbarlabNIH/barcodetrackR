options(shiny.maxRequestSize=1000*1024^2)

shinyServer(

  function(input, output) {
    output$rowsofmonks <- renderUI({
      dummyinteger <- as.numeric(input$nummonkeys)
      lapply(1:dummyinteger, function(i){
        fluidRow(
          column(12,
                 wellPanel(style="height:100px",
                   column(3,
                          textInput(paste0("monkey_name", i), label = paste0("Monkey Name ", i), value = "")
                   ),
                   column(3,
                          fileInput(paste0("monkey_out", i), label = paste0("Monkey Outfile ", i))

                   ),
                   column(3,
                          fileInput(paste0("monkey_key", i), label = paste0("Monkey Keyfile ", i))
                   ),
                   column(3,
                            numericInput("thresholdvalue", "Threshold", value = 2000)
                   )
                 )
          )
        )

      })
    })




  }
)



