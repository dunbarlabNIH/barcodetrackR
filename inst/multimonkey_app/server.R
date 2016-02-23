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
        tabPanel("Barcode Count",
                 uiOutput("barcodecount")),
        tabPanel("Diversity",
                 uiOutput("diversity")),
        tabPanel("Library Contributions",
                 uiOutput("lib_contrib")),
        tabPanel("Rank Abundance",
                 uiOutput("rankabundance"))
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

    data_as_list_subsets <- reactive({
      listy = list()
      for(i in 1:as.numeric(input$nummonkeys)){
        listy[[i]] = data_as_list()[[i]][,sample_list()[[i]],drop = FALSE]
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

    barcode_count_plotInput <- function(){
      print(barcodetrackR::multimonkey_richness_plot(
        outfile_list = data_as_list_subsets(),
        outfile_names = sample_names(),
        outfile_months_list = months_list(),
        celltypes_list = cells_list(),
        threshold_list = thresh_list(),
        combine = input$barcodecount_Combine,
        richness_type = input$barcodecount_Method,
        point_size = input$barcodecount_Dotsize,
        line_size = input$barcodecount_Linesize,
        y_lower = input$barcodecount_yLower,
        y_upper = input$barcodecount_yUpper))

    }

    output$barcode_count_plot <- renderPlot({
      barcode_count_plotInput()
      height = 500
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
                 selectInput("barcodecount_Method", label = "2. Barcode Plot Type", choices = c("cumulative", "unique", "new"), selected = "cumulative"),
                 selectInput("barcodecount_Combine", label = "3. Combine Monthly", choices = c(TRUE, FALSE)),
                 numericInput("barcodecount_Linesize", label = "4. Size of Line", value = 1.5, min = 0),
                 numericInput("barcodecount_Dotsize", label = "5. Size of Dot", value = 1.5, min = 0),
                 numericInput("barcodecount_yLower", label = "6. Lower Limit of Y axis", value = 0, min = 0),
                 numericInput("barcodecount_yUpper", label = "7. Upper Limit of Y axis", value = 3000, min = 0)
               )
        ),
        column(9,
               plotOutput('barcode_count_plot', height = 800))
      )
    })

    #========================================================================================================================
    #RANKABUNDANCE STUFF

    rankabundance_totaldata <- reactive({
      col_lengths <- unlist(lapply(data_as_list(), nrow))
      max_rows <- max(col_lengths)
      listy = list()
      for(i in 1:as.numeric(input$nummonkeys)){
        listy[[i]] = data_as_list()[[i]][1:max_rows,,drop = FALSE]
        listy[[i]][is.na(listy[[i]])] <- 0
      }
      temp_data <- do.call(cbind, listy)
      rownames(temp_data) <- NULL
      return(temp_data)
    })

    rankabundance_data <- reactive({
      return(rankabundance_totaldata()[,input$rankabundance_Samples,drop = FALSE])
    })

    rankabundance_plotInput <- function(){
      print(barcodetrackR::rank_abundance_plot(rankabundance_data(), dot_size = input$rankabundance_Dotsize, text_size = input$rankabundance_Textsize))
    }

    output$rankabundance_plot <- renderPlot({
      rankabundance_plotInput()
      height = 800
    })



    output$rankabundance <- renderUI({
      fluidRow(
        column(3,
               wellPanel(
                 selectInput("rankabundance_Samples", label = "1. Which Samples to Use", choices = colnames(rankabundance_totaldata()), multiple = TRUE),
                 numericInput("rankabundance_Dotsize", "2. Enter Dot Size: ", value = 3),
                 numericInput("rankabundance_Textsize", "3. Enter Text Size: ", value = 15)
               )),
        column(9,
               plotOutput("rankabundance_plot", height = 800))
      )

    })

  }
)



