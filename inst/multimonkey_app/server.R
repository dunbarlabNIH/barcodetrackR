options(shiny.maxRequestSize=1000*1024^2)

shinyServer(
  function(input, output) {
    #========================================================================================================================
    output$rowsofmonks <- renderUI({
      dummyinteger <- as.numeric(input$nummonkeys)
      lapply(1:dummyinteger, function(i){
        fluidRow(
          column(12,
                 wellPanel(style="height:100px",
                           column(3, textInput(paste0("monkey_name", i), label = paste0("Monkey Name ", i), value = "")),
                           column(3, fileInput(paste0("monkey_out", i), label = paste0("Monkey Outfile ", i))),
                           column(3, fileInput(paste0("monkey_key", i), label = paste0("Monkey Keyfile ", i))),
                           column(3, numericInput(paste0("monkey_threshold",i), "Threshold", value = 2000))
                 )
          )
        )
      })
    })
    #========================================================================================================================

    data_as_list <- eventReactive(input$start_button, {
      withProgress(message = "Applying Threshold", value = 0, {
        mylist <- list()
        for(i in c(1:as.numeric(input$nummonkeys))){
          incProgress(1/(as.numeric(input$nummonkeys)),paste0("Applying Thresh to Monkey ", i))
          mylist[[i]] <- barcodetrackR::threshold(read.delim(file = input[[paste0("monkey_out", i)]]$datapath, row.names =1),
                                                  input[[paste0("monkey_threshold",i)]])
          keyfile <- read.delim(input[[paste0("monkey_key", i)]]$datapath, row.names = 1)
          keyfile = keyfile[match(colnames(mylist[[i]]), rownames(keyfile)),]
          colnames(mylist[[i]]) <- paste(input[[paste0("monkey_name",i)]],keyfile$GIVENNAME)
        }
      })
      return(mylist)
    })
    #========================================================================================================================

    output$panel <- renderUI({
      if (is.null(data_as_list()))
        return()
      tabsetPanel(
        tabPanel("Barcode Count & Diversity",
                 uiOutput("barcodecount")),
        tabPanel("Library Contributions",
                 uiOutput("lib_contrib")),
        tabPanel("Rank Abundance",
                 uiOutput("rank_abundance"))
      )
    })

    #========================================================================================================================
    #BARCODECOUNT STUFF



    barcodecount_Modals <- function(){
      bM = lapply(1:(as.numeric(input$nummonkeys)),function(i){
        shinyBS::bsModal(paste0("barcodecount_",i), "Upload Data", paste0("barcodecount_Button",i), size = "large",
                         fluidRow(
                           column(12, wellPanel(
                             selectInput(paste0("barcodecount_Samples",i), label = "1. Which Samples to Use (order by month, then by celltype: e.g. 1m T, 1m B, 1m Gr, 2m T, 2m B, etc.)", choices = colnames(data_as_list()[[i]]), multiple = TRUE),
                             textInput(paste0("barcodecount_Months", i), '2. Enter months (in order) seperated by a comma: ', value = ""),
                             textInput(paste0("barcodecount_Celltypes",i), '3. Enter celltypes (in order) seperated by a comma: ', value = ""),
                             numericInput(paste0("barcodecount_Threshold",i), "4. Enter threshold: ", value = 0, min = 0, step = 1000))
                           )
                         )
        )
      })
    }

    barcodecount_Buttons <- function(){
      bc_BT =lapply(1:(as.numeric(input$nummonkeys)), function(i){
        actionButton(paste0("barcodecount_Button",i), paste0("Data set up for: ",input[[paste0("monkey_name", i)]]))
      })
      return(bc_BT)
    }

    sample_names <- reactive({
      listy = list()
      for(i in 1:as.numeric(input$nummonkeys)){
        listy[[i]] = input[[paste0("monkey_name",i)]]
      }
      return(listy)
    })

    sample_list <- reactive({
      listy = list()
      for(i in 1:as.numeric(input$nummonkeys)){
        listy[[i]] = input[[paste0("barcodecount_Samples",i)]]
      }
      return(listy)
    })

    months_list <- reactive({
      listy = list()
      for(i in 1:as.numeric(input$nummonkeys)){
        listy[[i]] = as.numeric(unlist(strsplit(input[[paste0("barcodecount_Months",i)]], split = ',')))
      }
      return(listy)
    })

    cells_list <- reactive({
      listy = list()
      for(i in 1:as.numeric(input$nummonkeys)){
        listy[[i]] = unlist(strsplit(input[[paste0("barcodecount_Celltypes", i)]], split = ','))
      }
      return(listy)
    })

    thresh_list <- reactive({
      listy = list()
      for(i in 1:as.numeric(input$nummonkeys)){
        listy[[i]] = as.numeric(input[[paste0("barcodecount_Threshold",i)]])
      }
      return(listy)
    })


    output$barcodecount <- renderUI({
      fluidRow(
        column(3,
               wellPanel(
                 strong("1. Data set up for monkeys: "),
                 br(),
                 barcodecount_Buttons(),
                 barcodecount_Modals(),
                 br(),
                 br(),
                 selectInput("barcodecount_Method", label = "2. Barcode Plot Type", choices = c("cumulative", "unique", "new", "diversity"), selected = "cumulative"),
                 selectInput("barcodecount_Combine", label = "3. Combine Monthly (if cumulative/unique/new chosen)", choices = c(TRUE, FALSE)),
                 selectInput("barcodecount_Diversitytype", label = "4. Diversity Type (if diversity chosen)", choices = c("shannon", "simpson", "invsimpson", selected = "shannon")),
                 numericInput("barcodecount_Linesize", label = "5. Size of Line", value = 1.5, min = 0),
                 numericInput("barcodecount_Dotsize", label = "6. Size of Dot", value = 1.5, min = 0),
                 numericInput("barcodecount_yLower", label = "7. Lower Limit of Y axis", value = 0, min = 0),
                 numericInput("barcodecount_yUpper", label = "8. Upper Limit of Y axis", value = 3000, min = 0)
               )
        ),
        column(9,
               plotOutput(plot(mtcars), height = 700)
        )
      )
    })








  }
)



