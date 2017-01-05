options(shiny.maxRequestSize=1000*1024^2)

shinyServer(

  function(input, output, session) {

    #=========================================================================================================

    output$thresholdPanel <- renderUI({
      if (is.null(input$file1) | is.null(input$file2) | is.null(input$file3))
        return()
      strong("5. Press Button to Apply Threshold")
      actionButton("threshybutton", "Apply Threshold",width="100%")


    })

    thresholded_data <- eventReactive(input$threshybutton, {
      withProgress(message = "Applying Threshold", value = 0, {
        df <- barcodetrackR::threshold(read.delim(input$file1$datapath, row.names = 1), input$thresholdvalue)
        incProgress(0.5, detail = "Applying Threshold")
        tdf <- t(df)
        keyfile <- read.delim(input$file2$datapath, stringsAsFactors = FALSE)
        keyfile <- keyfile[,c("FILENAME", "GIVENNAME")]
        if(any(duplicated(keyfile$FILENAME))){
          stop("Keyfile contains duplicate FILENAME")
        }
        if(any(duplicated(keyfile$GIVENNAME))){
          stop("Keyfile contains duplicate GIVENNAME")
        }
        if(!(all(rownames(tdf) %in% keyfile$FILENAME))){
          stop("Missing certain files in keyfile")
        }
        if(!(all(keyfile$FILENAME %in% rownames(tdf)))){
          stop("FILENAME in keyfile is not in outfile")
        }
        if(length(setdiff(keyfile$FILENAME, rownames(tdf))) != 0){
          print(setdiff(keyfile$FILENAME, rownames(tdf)))
          stop("Number of samples in keyfile differs from number of samples in outfile")
        }
        keyfile <- keyfile[match(rownames(tdf), keyfile$FILENAME),]
        df <- data.frame(GIVENNAME = keyfile$GIVENNAME, tdf)
      })
      df <- data.frame(df)
      return(df)
    }
    )


    current_threshold <- eventReactive(input$threshybutton,{
      return(paste("Current threshold applied is:\n", input$thresholdvalue, sep = ""))
    })


    output$thresholdInfo <- renderUI({
      if (is.null(thresholded_data()))
        return()
      strong(current_threshold())
    })

    readme_data <- eventReactive(input$threshybutton,{
      readme <- read.delim(input$file3$datapath, skip = 4, header = FALSE, stringsAsFactors = FALSE)
      readme <- readme[1:(nrow(readme)-1),1:3]
      colnames(readme) <- c("FILENAME", "MAPPED", "READS")
      keyfile <- read.delim(input$file2$datapath, stringsAsFactors = FALSE)
      readme$PERCENTAGE <- readme$MAPPED/readme$READS * 100
      readme$GIVENNAME <- keyfile$GIVENNAME[match(readme$FILENAME, keyfile$FILENAME)]
      return(readme)
    })



    #======================================================================================================

    #TABPANEL

    output$Panel <- renderUI({

      if (is.null(thresholded_data()))
        return()

      tabsetPanel(
        tabPanel("DataStatistics",
                 uiOutput("DataStatistics")),
        tabPanel("Heatmap",
                 uiOutput("Heatmap")),
        tabPanel("CorPlot",
                 uiOutput("CorPlot")),
        tabPanel("Diversity & Richness",
                 uiOutput("Diversity")),
        tabPanel("Scatter Plot",
                 uiOutput("ScatterPlot")),
        tabPanel("Tree Map",
                 uiOutput("TreeMap")),
        tabPanel("BBHM",
                 uiOutput("BBHM")),
        tabPanel("TernPlot",
                 uiOutput("TernPlot")),
        tabPanel("RankAbundance",
                 uiOutput("RankAbundance")),
        tabPanel("UnilineageBias",
                 uiOutput("UnilineageBias")),
        tabPanel("ClonalBias",
                 uiOutput("ClonalBias"))
      )
    })


    #======================================================================================================

    #DATASTATISTICS TAB

    output$DataStatistics <- renderUI({
      fluidRow(
        column(7,dataTableOutput('renderedReadme')),
        column(5, plotOutput('readmeHistogram'))
      )
    })

    output$readmeHistogram <- renderPlot({
      hist(readme_data()$MAPPED/readme_data()$READS * 100,
           breaks = seq(0,100, by = 5),
           xlim = c(0,100),
           main = "MAPPING % HISTOGRAM",
           col = "lightblue",
           xlab = "PERCENTAGE MAPPED")
    })
    output$renderedReadme <- renderDataTable(readme_data(), options = list(scrollX=TRUE, pageLength = 10))


    #======================================================================================================

    #HEATMAP TAB

    output$Heatmap <- renderUI({

      if (is.null(thresholded_data()))
        return()

      HeatmapInput <- function(){
        print(barcodetrackR::barcode_ggheatmap(your_data = Heatmap_data(),
                                               n_clones = input$Heatmap_top_clones,
                                               your_title = input$Heatmap_title,
                                               grid = input$Heatmap_grid,
                                               label_size = input$Heatmap_labels,
                                               dendro = input$Heatmap_dendrogram,
                                               cellnote_size = input$Heatmap_starsize,
                                               log_transform = input$Heatmap_log_transform,
                                               log_choice = switch(as.character(input$Heatmap_scale), "2" = 2, "e" = exp(1), "10" = 10, "100" = 100),
                                               distance_method = input$Heatmap_distance,
                                               minkowski_power = input$Heatmap_mink_distance,
                                               cellnote_option = input$Heatmap_cellnote_option,
                                               hclust_linkage = input$Heatmap_hclust_linkage,
                                               row_order = input$Heatmap_row_order,
                                               clusters = input$Heatmap_clusters))}

      HeatmapREADMEtable <- function(){
        your_table <- readme_data()
        your_table <- your_table[match(input$Heatmap_samples, your_table$GIVENNAME),]
        temp_names <- your_table$GIVENNAME
        temp_clones <- colSums(Heatmap_data() > 0)
        your_table <- as.data.frame(cbind(your_table$MAPPED, your_table$READS, round(your_table$MAPPED/your_table$READS * 100,2), temp_clones))
        rownames(your_table) <- temp_names
        colnames(your_table) <- c("MAPPED", "READS", "MAP %", "CLONES")
        return(cbind(rownames(t(your_table)),t(your_table)))
      }


      output$heatmap_datatable <- renderDataTable({HeatmapREADMEtable()},
                                                  options = list(paging = FALSE, searching = FALSE))


      output$downloadHeatmapkey <- downloadHandler(
        filename = function() {paste(input$file1, "_heatmapkey.txt", sep = "")},
        content = function(file){
          write.table(barcodetrackR::barcode_ggheatmap(your_data = Heatmap_data(),
                                                       n_clones = input$Heatmap_top_clones,
                                                       log_transform = input$Heatmap_log_transform,
                                                       printtable = TRUE,
                                                       table_option = input$Heatmap_table_option,
                                                       log_choice = switch(as.character(input$Heatmap_scale), "2" = 2, "e" = exp(1), "10" = 10),
                                                       distance_method = input$Heatmap_distance,
                                                       minkowski_power = input$Heatmap_mink_distance,
                                                       hclust_linkage = input$Heatmap_hclust_linkage), file, sep = '\t', quote = FALSE)
        }
      )

      output$downloadHeatmapSTATS <- downloadHandler(
        filename = function() {paste(input$file1, "_heatmapSTATS.txt", sep = "")},
        content = function(file){
          write.table(HeatmapREADMEtable(), file, sep = '\t', quote = FALSE)
        }
      )

      output$viewHeatmap <- renderPlot({
        HeatmapInput()
      })


      Heatmap_data <- reactive({
        df <- thresholded_data()
        df <- df[df$GIVENNAME %in% input$Heatmap_samples,] #subset samples
        df$GIVENNAME <- factor(df$GIVENNAME, levels = input$Heatmap_samples)
        df <- df[order(df$GIVENNAME),]
        newcolnames <- df$GIVENNAME
        df$GIVENNAME <- NULL
        df <- data.frame(t(df))
        colnames(df) <- newcolnames
        return(df)

      })

      output$Heatmap_download_samples <- downloadHandler(
        filename = function() {paste("_samplelist.txt", sep = "\t")},
        content = function(file){
          write.table(input$Heatmap_samples, file, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
        }
      )

      Heatmap_uploaded_samples <- reactive({
        samples <- as.vector(t(read.delim(input$Heatmap_uploaded_samples$datapath, stringsAsFactors = FALSE)))
        return(samples)
      })

      observeEvent(input$Heatmap_uploaded_samples, {updateSelectizeInput(session, inputId = 'Heatmap_samples', selected = Heatmap_uploaded_samples())})

      fluidRow(
        column(3,
               wellPanel(
                 div(style="display:inline-block; height:85px;",fileInput("Heatmap_uploaded_samples", "1. Upload Prepared Sample List or Input Samples")),
                 selectizeInput("Heatmap_samples", label = NULL, choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = TRUE),
                 downloadButton("Heatmap_download_samples", "Download This Sample List"),
                 br(),
                 br(),
                 numericInput("Heatmap_top_clones", "2. Number of top clones", value = 10),
                 strong("3. Options"),
                 checkboxInput("Heatmap_grid", label = "Grid", value = TRUE),
                 checkboxInput("Heatmap_log_transform", label = "Log-Transform", value = TRUE),
                 checkboxInput("Heatmap_dendrogram", label = "Display Dendrogram", value = FALSE),
                 selectInput("Heatmap_cellnote_option", "4. Select Cell Display Option",
                             choices = c("reads", "percents", "logs", "stars", "ranks"),
                             selected = "stars"),
                 numericInput("Heatmap_labels", "5. Set Column Label Size", value = 20),
                 numericInput("Heatmap_starsize", "6. Set Cell Label Size", value = 15),
                 textInput("Heatmap_title", "7. Title for Heatmap", value = ""),
                 selectInput("Heatmap_distance", "8. Select Distance Metric/Function",
                             choices = sort(as.vector(unlist(summary(proxy::pr_DB)[1]))),
                             selected = "Euclidean"),
                 numericInput("Heatmap_mink_distance", "9. If Minkowski, choose Minkowski Power", value = 2, step = 1),
                 selectInput("Heatmap_hclust_linkage", "10. Select Clustering Linkage",
                             choices = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"),
                             selected = "complete"),
                 selectInput("Heatmap_scale", "11. Select Log",
                             choices = c("2", "e", "10", "100"),
                             selected = "e"),
                 numericInput("Heatmap_clusters", "12. Select number of clusters to cut", value = 0, step = 1, min = 1, max = 11),
                 selectInput("Heatmap_row_order", "13. How to order rows", choices = c("hierarchical", "emergence"), selected = "hierarchical"),
                 selectInput("Heatmap_table_option", "14. Select Format for Key (below)",
                             choices = c("logs", "reads", "percents", "ranks"),
                             selected = "percents"),
                 strong("15. Press button to download Heatmap Key."),
                 br(),
                 downloadButton('downloadHeatmapkey', 'Heatmap_key'),
                 br(),
                 strong("16. Press button to download Heatmap Stats."),
                 br(),
                 downloadButton('downloadHeatmapSTATS', 'Heatmap_STATS')
               )
        ),
        column(8,
               plotOutput('viewHeatmap', height = 800),
               dataTableOutput('heatmap_datatable')
        )
      )
    })

    #======================================================================================================

    #CORPLOT TAB

    output$CorPlot <- renderUI({


      if (is.null(thresholded_data()))
        return()

      observeEvent(input$corplot_Copier, {updateSelectizeInput(session, inputId = 'corplot_Samples', selected = input$Heatmap_samples)})


      corplotInput <- function(){
        barcodetrackR::cor_plot(your_data = corplot_data(),
                                names = colnames(corplot_data()),
                                thresh = input$corplot_thresh,
                                your_title = input$corplot_Title,
                                method_corr = input$corplot_Method,
                                labelsizes = input$corplot_Labels,
                                plottype = input$corplot_Type,
                                no_negatives = input$corplot_excludeneg,
                                show_grid = input$corplot_Grid,
                                colorscale = input$corplot_Colors)
      }

      output$downloadcorplotzip <- downloadHandler(
        filename = function() {"corplot_files.zip"},
        content = function(file){
          tmpdir <- tempdir()
          setwd(tmpdir)
          listoffiles <- barcodetrackR::cor_plot(your_data = corplot_data(), names = colnames(corplot_data()),
                                                 thresh = input$corplot_thresh, your_title = input$corplot_Title,
                                                 method_corr = input$corplot_Method, labelsizes = input$corplot_Labels,
                                                 plottype = input$corplot_Type, printtables = TRUE)
          for(i in seq_along(listoffiles)){
            write.table(listoffiles[[i]], file = paste0(names(listoffiles)[i], ".txt"), quote = FALSE, sep = '\t')
          }
          zip(zipfile = file, files = c("cortable.txt","cortable_pval.txt", "cortable_ci_hi.txt","cortable_ci_lo.txt"))
          if(file.exists(file)) {file.rename(paste0(file), file)}
        },
        contentType = "application/zip"
      )



      output$viewcorplot <- renderPlot({
        corplotInput()
      })

      corplot_data <- reactive({
        cf <- thresholded_data()
        cf <- cf[cf$GIVENNAME %in% input$corplot_Samples,] #subset samples
        cf$GIVENNAME <- factor(cf$GIVENNAME, levels = input$corplot_Samples)
        cf <- cf[order(cf$GIVENNAME),]
        newcolnames <- cf$GIVENNAME
        cf$GIVENNAME <- NULL
        cf <- data.frame(t(cf))
        colnames(cf) <- newcolnames
        return(cf)

      })

      #======================================================================================================


      fluidRow(
        column(3,
               wellPanel(
                 selectizeInput("corplot_Samples", label = "1. Which Samples to Use (order matters)",
                                choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = TRUE, selected = ""),
                 actionButton("corplot_Copier", label = "Copy BCHM Samples"),
                 br(),
                 br(),
                 numericInput("corplot_thresh", "2. Reads threshold", value = 0),
                 textInput("corplot_Title", "3. Title", value = ""),
                 selectInput("corplot_Type", '4. Choose Type of Plot', choices = c("circle", "square", "ellipse", "color", "number", "shade", "pie"), selected = "square"),
                 strong("5. Options"),
                 checkboxInput("corplot_excludeneg", "Exclude negatives", value = FALSE),
                 checkboxInput("corplot_Grid", "Grid ON/OFF", value = TRUE),
                 selectInput("corplot_Method", "6. Chooose Correlation Method", choices = c("pearson", "kendall", "spearman", "manhattan"), selected = "pearson"),
                 selectInput("corplot_Colors", "7. Choose Color Scale", choices = c("default", "rainbow", "white_heat"), selected = "default"),
                 numericInput("corplot_Labels", "8. Set Label Size", value = 2),
                 strong("9. Press button to downlaod CorPlot Files as .zip"),
                 br(),
                 downloadButton('downloadcorplotzip', 'corr_plot.zip')
               )
        ),

        column(8,
               plotOutput('viewcorplot', height = 800)
        )
      )
    })


    #======================================================================================================
    #TOPCLONESTRACKER TAB

    output$TopClonesTracker <- renderUI({
      if(is.null(thresholded_data()))
        return()

      topclonestrackerInput <- function(){
        print(barcodetrackR::topclonestracker(topclonestracker_data(),
                                              top_clones_choice = topclonestracker_Choice(),
                                              months = topclonestracker_Months(),
                                              n_clones = input$topclonestracker_Clones,
                                              linesize = input$topclonestracker_Linesize,
                                              text_size = input$topclonestracker_Textsize,
                                              y_limit = input$topclonestracker_yLim,
                                              plot_theme = input$topclonestracker_Theme,
                                              your_title = input$topclonestracker_Title))

      }

      output$viewtopclonestracker <- renderPlot({
        topclonestrackerInput()
        height = 700
      })

      topclonestracker_data <- reactive({
        ff <- thresholded_data()
        ff <- ff[ff$GIVENNAME %in% input$topclonestracker_Samples,] #subset samples
        ff$GIVENNAME <- factor(ff$GIVENNAME, levels = input$topclonestracker_Samples)
        ff <- ff[order(ff$GIVENNAME),]
        newcolnames <- ff$GIVENNAME
        ff$GIVENNAME <- NULL
        ff$EXPERIMENT <- NULL
        ff$CELLTYPE <- NULL
        ff$MONTH <- NULL
        ff$LOCATION <- NULL
        ff$MISC <- NULL
        ff <- data.frame(t(ff))
        colnames(ff) <- newcolnames
        return(ff)
      })

      topclonestracker_Choice <- reactive({
        tcchoice <- thresholded_data()
        tcchoice <- tcchoice[tcchoice$GIVENNAME == input$topclonestracker_Selected,,drop = FALSE]
        newcolnames <- tcchoice$GIVENNAME
        tcchoice$GIVENNAME <- NULL
        tcchoice$EXPERIMENT <- NULL
        tcchoice$CELLTYPE <- NULL
        tcchoice$MONTH <- NULL
        tcchoice$LOCATION <- NULL
        tcchoice$MISC <- NULL
        tcchoice <- data.frame(t(tcchoice))
        colnames(tcchoice) <- newcolnames
        return(tcchoice)
      })



      topclonestracker_Months <- reactive({
        return(as.numeric(unlist(strsplit(input$topclonestracker_Months, split = ','))))
      })

      fluidRow(
        column(3,
               wellPanel(
                 selectInput("topclonestracker_Samples", label = "1. Which Samples to Use",
                             choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = TRUE),
                 selectInput("topclonestracker_Selected", label = "2. Which Sample to Use for Top Clones",
                             choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = FALSE),
                 textInput("topclonestracker_Months", '3. Enter months (in order) seperated by a comma: ', value = ""),
                 numericInput("topclonestracker_Clones", "4. Select Number of Clones", value = 10),
                 numericInput("topclonestracker_Linesize", "5. Enter Line Size: ", value = 1),
                 numericInput("topclonestracker_Textsize", "6. Enter Text Size: ", value = 15),
                 numericInput("topclonestracker_yLim", "7. Enter Y Limit: ", value = 100),
                 selectInput("topclonestracker_Theme", "8. Choose Theme",
                             choices = c("original", "BW", "classic"), selected = "classic"),
                 textInput("topclonestracker_Title", "9. Title", value = "Stacked Area Plot")
               )),



        column(9,
               plotOutput('viewtopclonestracker', height = 700)
        )

      )






    })

    #======================================================================================================
    #DIVERSITY TAB

    output$Diversity <- renderUI({
      if(is.null(thresholded_data()))
        return()

      diversityInput <- function(){
        if (input$diversity_Type == "richness"){
          print(barcodetrackR::richness_plot(your_data = diversity_data(),
                                             months = diversity_Months(), celltypes = diversity_Celltypes(),
                                             thresh = input$diversity_Thresh, point_size = input$diversity_Pointsize,
                                             line_size = input$diversity_Linesize, richness_type = input$richness_Type,
                                             y_lower = input$diversity_yLower, y_upper = input$diversity_yUpper))

        } else {
          print(barcodetrackR::diversity_plot(your_data = diversity_data(),
                                              months = diversity_Months(), celltypes = diversity_Celltypes(),
                                              thresh = input$diversity_Thresh, point_size = input$diversity_Pointsize,
                                              line_size = input$diversity_Linesize, index_type = input$diversity_Type,
                                              y_lower = input$diversity_yLower, y_upper = input$diversity_yUpper))
        }
      }

      output$viewdiversity <- renderPlot({
        diversityInput()
        height = 700
      })

      diversity_data <- reactive({
        gf <- thresholded_data()

        gf <- gf[gf$GIVENNAME %in% input$diversity_Samples,] #subset samples
        gf$GIVENNAME <- factor(gf$GIVENNAME, levels = input$diversity_Samples)
        gf <- gf[order(gf$GIVENNAME),]
        newcolnames <- gf$GIVENNAME

        gf$GIVENNAME <- NULL
        gf$EXPERIMENT <- NULL
        gf$CELLTYPE <- NULL
        gf$MONTH <- NULL
        gf$LOCATION <- NULL
        gf$MISC <- NULL

        gf <- data.frame(t(gf))
        colnames(gf) <- newcolnames
        return(gf)

      })

      diversity_Months <- reactive({
        return(as.numeric(unlist(strsplit(input$diversity_Months, split = ','))))
      })

      diversity_Celltypes <- reactive({
        return(as.character(unlist(strsplit(input$diversity_Celltypes, split = ','))))
      })

      fluidRow(
        column(3,
               wellPanel(
                 selectInput("diversity_Samples", label = "1. Which Samples to Use (order by month, then by celltype: e.g. 1m T, 1m B, 1m Gr, 2m T, 2m B, etc.)",
                             choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = TRUE),

                 textInput("diversity_Months", '2. Enter months (in order) seperated by a comma: ', value = ""),

                 textInput("diversity_Celltypes", '3. Enter celltypes (in order) seperated by a comma: ', value = ""),

                 numericInput("diversity_Thresh", "4. Select Threshold", value = 0, step = 100, min = 0),

                 selectInput("diversity_Type", label  = "5A. Select Diversity Index", choices = c("shannon", "simpson", "invsimpson", "herfindahl-hirschman", "richness"), selected = "shannon"),

                 selectInput("richness_Type", label = "5B. If \"richness\" selected, please select type:", choices = c("unique", "cumulative", "new")),

                 numericInput("diversity_Pointsize", "6. Point Size: ", value = 3),

                 numericInput("diversity_Linesize", "7. Line Size: ", value = 2),

                 numericInput("diversity_yUpper", "8. y-upper Lim", value = 10),

                 numericInput("diversity_yLower", "9. y-lower Lim", value = 0)

               )),



        column(9,
               plotOutput('viewdiversity', height = 700)
        )

      )






    })





    #======================================================================================================
    #SCATTER PLOT TAB

    output$ScatterPlot <- renderUI({
      if(is.null(thresholded_data()))
        return()

      scatterInput <- function(){
        normdata <- 100*prop.table(as.matrix(scatter_data()), margin = 2)
        plot(normdata, col = 'black', ylim = c(0,input$scatter_Ylim), xlim = c(0,input$scatter_Xlim), cex = input$scatter_Dotsize)
        legend("topright", legend = paste("Correlation (R) =", cor(normdata[,1], normdata[,2])))
      }

      output$viewscatter <- renderPlot({
        scatterInput()
        height = 700
      })

      scatter_data <- reactive({
        scf <- thresholded_data()

        scf <- scf[scf$GIVENNAME %in% c(input$scatter_Sample1, input$scatter_Sample2),] #subset samples
        newcolnames <- scf$GIVENNAME

        scf$GIVENNAME <- NULL
        scf$EXPERIMENT <- NULL
        scf$CELLTYPE <- NULL
        scf$MONTH <- NULL
        scf$LOCATION <- NULL
        scf$MISC <- NULL

        scf <- data.frame(t(scf))
        colnames(scf) <- newcolnames
        return(scf)

      })

      fluidRow(
        column(3,
               wellPanel(
                 selectInput("scatter_Sample1", label = "1. First Sample to Use", choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = FALSE),
                 selectInput("scatter_Sample2", label = "2. Second Sample to Use", choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = FALSE),
                 numericInput("scatter_Threshold", "3. Select Threshold", value = 0, step = 100, min = 0),
                 numericInput("scatter_Ylim", "4. Y Limit: ", value = 100, step = 0.1, min = 0),
                 numericInput("scatter_Xlim", "5. X Limit: ", value = 100, step = 0.1, min = 0),
                 numericInput("scatter_Dotsize","6. Dot Size: ", value = 1)

               )),



        column(7,
               plotOutput('viewscatter', height = 700)
        )

      )


    })

    #======================================================================================================
    #TREEMAP TAB

    output$TreeMap <- renderUI({

      if(is.null(thresholded_data()))
        return()

      output$viewtreemap <- renderPlot({
        barcodetrackR::barcode_treemap(your_data = treemap_data(), column_choice = colnames(treemap_data()),
                                       threshold = input$treemap_Threshold,
                                       colors_choice = c(input$treemap_LoColor, input$treemap_HiColor),
                                       your_title = input$treemap_Title)
        height = 700
      })

      treemap_data <- reactive({
        tmdf <- thresholded_data()
        tmdf <- tmdf[tmdf$GIVENNAME %in% input$treemap_Sample,,drop = FALSE] #subset samples
        newcolnames <- tmdf$GIVENNAME
        tmdf$GIVENNAME <- NULL
        tmdf$EXPERIMENT <- NULL
        tmdf$CELLTYPE <- NULL
        tmdf$MONTH <- NULL
        tmdf$LOCATION <- NULL
        tmdf$MISC <- NULL
        tmdf <- data.frame(t(tmdf))
        colnames(tmdf) <- newcolnames
        return(tmdf)
      })

      fluidRow(
        column(3,
               wellPanel(
                 selectInput("treemap_Sample", label = "1. Which Sample to Use",
                             choices = as.vector(unique(thresholded_data()$GIVENNAME))),

                 numericInput("treemap_Threshold", "2. Select Threshold", value = 0, step = 100, min = 0),

                 selectInput("treemap_HiColor", label  = "3. Select Upper Color",
                             choices = c("red", "orange", "yellow", "green", "blue", "violet", "black", "white", "grey"),
                             selected = "yellow"),

                 selectInput("treemap_LoColor", label = "4. Select Lower Color",
                             choices = c("red", "orange", "yellow", "green", "blue", "violet", "black", "white", "grey"),
                             selected = "blue"),
                 textInput("treemap_Title", "5. Title", value = "")
               )),
        column(9,
               plotOutput('viewtreemap', height = 700)
        )

      )



    })

    #======================================================================================================
    #BBHM TAB

    output$BBHM <- renderUI({


      if (is.null(thresholded_data()))
        return()


      BBHMInput <- function(){
        barcodetrackR::BBHM(your_data = BBHM_data(),
                            col_labels = input$BBHM_labels,
                            threshold = input$BBHM_threshold
        )
      }

      output$viewBBHM <- renderPlot({
        BBHMInput()
        height = 900
      })


      BBHM_data <- reactive({
        BBdf <- thresholded_data()
        BBdf <- BBdf[BBdf$GIVENNAME %in% input$BBHM_samples,] #subset samples
        BBdf$GIVENNAME <- factor(BBdf$GIVENNAME, levels = input$BBHM_samples)
        BBdf <- BBdf[order(BBdf$GIVENNAME),]
        newcolnames <- BBdf$GIVENNAME
        BBdf$GIVENNAME <- NULL
        BBdf$EXPERIMENT <- NULL
        BBdf$CELLTYPE <- NULL
        BBdf$MONTH <- NULL
        BBdf$LOCATION <- NULL
        BBdf$MISC <- NULL
        BBdf <- data.frame(t(BBdf))
        colnames(BBdf) <- newcolnames
        return(BBdf)
      })

      fluidRow(
        column(3,
               wellPanel(
                 selectizeInput("BBHM_samples", label = "1. Which Samples to Use (in order)",
                                choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = TRUE),
                 numericInput("BBHM_threshold", "2. Set Threshold", value = 0),
                 numericInput("BBHM_labels", "3. Set Column Label Size", value = 1.5)
               )
        ),
        column(9,
               plotOutput('viewBBHM', height = 700)
        )

      )

    })

    #======================================================================================================
    #TERNPLOT TAB

    output$TernPlot <- renderUI({
      if(is.null(thresholded_data()))
        return()

      ternplotInput <- function(){
        print(barcodetrackR::ternary_plot(ternplot_data(),
                                          dot_size = input$ternplot_Dotsize,
                                          show_arrows = input$ternplot_Showarrows,
                                          show_ticks = input$ternplot_Showticks,
                                          density_mode = input$ternplot_Density))

      }

      output$viewternplot<- renderPlot({
        ternplotInput()
        height = 700
      })

      ternplot_data <- reactive({
        terndf <- thresholded_data()

        terndf <- terndf[terndf$GIVENNAME %in% input$ternplot_Samples,] #subset samples
        terndf$GIVENNAME <- factor(terndf$GIVENNAME, levels = input$ternplot_Samples)
        terndf <- terndf[order(terndf$GIVENNAME),]
        newcolnames <- terndf$GIVENNAME

        terndf$GIVENNAME <- NULL
        terndf$EXPERIMENT <- NULL
        terndf$CELLTYPE <- NULL
        terndf$MONTH <- NULL
        terndf$LOCATION <- NULL
        terndf$MISC <- NULL

        terndf <- data.frame(t(terndf))
        colnames(terndf) <- newcolnames
        return(terndf)

      })

      fluidRow(
        column(3,
               wellPanel(
                 selectInput("ternplot_Samples", label = "1. Which Samples to Use (3)",
                             choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = TRUE),
                 strong("2. Options"),
                 checkboxInput("ternplot_Density", label = "Density Mode", value = FALSE),
                 checkboxInput("ternplot_Showarrows", label = "Show Arrows", value = TRUE),
                 checkboxInput("ternplot_Showticks", label = "Show Tick Marks", value = FALSE),
                 numericInput("ternplot_Dotsize", "3. Enter Dot Size: ", value = 1000)
               )),



        column(9,
               plotOutput('viewternplot', height = 700)
        )

      )






    })


    #======================================================================================================
    #RANKABUNDANCE TAB

    output$RankAbundance <- renderUI({
      if(is.null(thresholded_data()))
        return()

      rankabundanceInput <- function(){
        print(barcodetrackR::rank_abundance_plot(rankabundance_data(),
                                                 dot_size = input$rankabundance_Dotsize,
                                                 text_size = input$rankabundance_Textsize))

      }

      output$viewrankabundance<- renderPlot({
        rankabundanceInput()
        height = 700
      })

      rankabundance_data <- reactive({
        rankabdf <- thresholded_data()
        rankabdf <- rankabdf[rankabdf$GIVENNAME %in% input$rankabundance_Samples,] #subset samples
        rankabdf$GIVENNAME <- factor(rankabdf$GIVENNAME, levels = input$rankabundance_Samples)
        rankabdf <- rankabdf[order(rankabdf$GIVENNAME),]
        newcolnames <- rankabdf$GIVENNAME

        rankabdf$GIVENNAME <- NULL
        rankabdf$EXPERIMENT <- NULL
        rankabdf$CELLTYPE <- NULL
        rankabdf$MONTH <- NULL
        rankabdf$LOCATION <- NULL
        rankabdf$MISC <- NULL

        rankabdf <- data.frame(t(rankabdf))
        colnames(rankabdf) <- newcolnames
        return(rankabdf)

      })

      fluidRow(
        column(3,
               wellPanel(
                 selectInput("rankabundance_Samples", label = "1. Which Samples to Use",
                             choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = TRUE),
                 numericInput("rankabundance_Dotsize", "2. Enter Dot Size: ", value = 3),
                 numericInput("rankabundance_Textsize", "3. Enter Text Size: ", value = 15)
               )),



        column(9,
               plotOutput('viewrankabundance', height = 700)
        )

      )






    })



    #======================================================================================================
    #UNILINEAGEBIAS TAB

    output$UnilineageBias <- renderUI({
      if(is.null(thresholded_data()))
        return()

      unilineagebiasInput <- function(){
        if(input$unilineagebias_plot_mode == "heatmap"){
          barcodetrackR::unilineage_bias(your_data = unilineagebias_data(),
                                         CT = input$unilineagebias_CT,
                                         TP = input$unilineagebias_TP,
                                         percent_thresh = input$unilineagebias_thresh,
                                         ratio_thresh = input$unilineagebias_ratio,
                                         line_months = unilineagebias_line_months(),
                                         only_biased = input$unilineagebias_only_biased,
                                         plot_mode = input$unilineagebias_plot_mode,
                                         text_size = input$unilineagebias_text_size,
                                         line_size = input$unilineagebias_line_size,
                                         dot_size = input$unilineagebias_dot_size,
                                         your_title = input$unilineagebias_your_title,
                                         y_upper = input$unilineagebias_y_upper,
                                         y_lower = input$unilineagebias_y_lower,
                                         cellnote_display = input$unilineagebias_cellnote_display,
                                         by_celltype = input$unilineagebias_by_celltype,
                                         cellnote_size = input$unilineagebias_cellnote_size,
                                         print_table = FALSE)

        } else {
          print(barcodetrackR::unilineage_bias(your_data = unilineagebias_data(),
                                               CT = input$unilineagebias_CT,
                                               TP = input$unilineagebias_TP,
                                               percent_thresh = input$unilineagebias_thresh,
                                               ratio_thresh = input$unilineagebias_ratio,
                                               line_months = unilineagebias_line_months(),
                                               only_biased = input$unilineagebias_only_biased,
                                               plot_mode = input$unilineagebias_plot_mode,
                                               text_size = input$unilineagebias_text_size,
                                               line_size = input$unilineagebias_line_size,
                                               dot_size = input$unilineagebias_dot_size,
                                               your_title = input$unilineagebias_your_title,
                                               y_upper = input$unilineagebias_y_upper,
                                               y_lower = input$unilineagebias_y_lower,
                                               cellnote_display = input$unilineagebias_cellnote_display,
                                               by_celltype = input$unilineagebias_by_celltype,
                                               cellnote_size = input$unilineagebias_cellnote_size,
                                               print_table = FALSE))
        }
      }

      output$viewunilineagebias<- renderPlot({
        unilineagebiasInput()
        height = 700
      })

      unilineagebias_data <- reactive({
        uni_df <- thresholded_data()
        uni_df <- uni_df[uni_df$GIVENNAME %in% input$unilineagebias_samples,] #subset samples
        uni_df$GIVENNAME <- factor(uni_df$GIVENNAME, levels = input$unilineagebias_samples)
        uni_df <- uni_df[order(uni_df$GIVENNAME),]
        newcolnames <- uni_df$GIVENNAME
        uni_df$GIVENNAME <- NULL
        uni_df$EXPERIMENT <- NULL
        uni_df$CELLTYPE <- NULL
        uni_df$MONTH <- NULL
        uni_df$LOCATION <- NULL
        uni_df$MISC <- NULL
        uni_df <- data.frame(t(uni_df))
        colnames(uni_df) <- newcolnames
        return(uni_df)
      })

      unilineagebias_line_months <- reactive({
        return(as.numeric(unlist(strsplit(input$unilineagebias_line_months, split = ','))))
      })

      fluidRow(
        column(3,
               wellPanel(
                 selectInput("unilineagebias_samples", label = "1. Which Samples to Use",
                             choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = TRUE),
                 numericInput("unilineagebias_CT", label = "2. Number of Cell Types", value = 5),
                 numericInput("unilineagebias_TP", label = "3. Number of Time Points", value = 1),
                 numericInput("unilineagebias_thresh", label = "4. Threshold of Bias", value = 0.01),
                 numericInput("unilineagebias_ratio", label = "5. Ratio or Bias", value = 10),
                 textInput("unilineagebias_line_months", '6. Enter months (in order) seperated by a comma: ', value = ""),
                 strong("7. Options"),
                 checkboxInput("unilineagebias_only_biased", label = "Only Biased", value = FALSE),
                 checkboxInput("unilineagebias_by_celltype", label = "By Celltype", value = FALSE),
                 selectInput("unilineagebias_plot_mode", label = "8. Choose plot type",
                             choices = c("heatmap", "bar", "line"), multiple = FALSE),
                 numericInput("unilineagebias_text_size", label = "9. Text Size", value = 10),
                 numericInput("unilineagebias_line_size", label = "10. Line Size", value = 3),
                 numericInput("unilineagebias_dot_size", label = "11. Dot Size", value = 5),
                 textInput("unilineagebias_your_title", label = "12. Your title", value = ""),
                 numericInput("unilineagebias_y_upper", label = "13. Upper Y limit", value = 1),
                 numericInput("unilineagebias_y_lower", label = "14. Lower Y Limit", value = 0),
                 selectInput("unilineagebias_cellnote_display", label = "15. Cellnote Labels",
                             choices = c("stars", "percents")),
                 numericInput("unilineagebias_cellnote_size", label = "16. Size of Cell Labels", value = 3)
               )),
        column(9,
               plotOutput('viewunilineagebias', height = 900)
        )
      )
    })

    #======================================================================================================
    #CLONALBIAS TAB

    output$ClonalBias <- renderUI({
      if(is.null(thresholded_data()))
        return()

      clonalbiasInput <- function(){
        if(input$clonalbias_Type == "Dots"){
          print(barcodetrackR::dot_bias(clonalbias_data(), text_size = input$clonalbias_Textsize))
        } else if(input$clonalbias_Type == "Bars"){
          print(barcodetrackR::clonal_bias(clonalbias_data(), text_size = input$clonalbias_Textsize))
        }
      }

      output$viewclonalbias<- renderPlot({
        clonalbiasInput()
        height = 700
      })

      clonalbias_data <- reactive({
        cb_df <- thresholded_data()
        cb_df <- cb_df[cb_df$GIVENNAME %in% input$clonalbias_Samples,] #subset samples
        cb_df$GIVENNAME <- factor(cb_df$GIVENNAME, levels = input$clonalbias_Samples)
        cb_df <- cb_df[order(cb_df$GIVENNAME),]
        newcolnames <- cb_df$GIVENNAME
        cb_df$GIVENNAME <- NULL
        cb_df$EXPERIMENT <- NULL
        cb_df$CELLTYPE <- NULL
        cb_df$MONTH <- NULL
        cb_df$LOCATION <- NULL
        cb_df$MISC <- NULL
        cb_df <- data.frame(t(cb_df))
        colnames(cb_df) <- newcolnames
        return(cb_df)
      })
      fluidRow(
        column(3,
               wellPanel(
                 selectInput("clonalbias_Samples", label = "1. Which Samples to Use",
                             choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = TRUE),
                 numericInput("clonalbias_Textsize", "2. Enter Text Size: ", value = 15),
                 selectInput("clonalbias_Type", label = "3. Pick type of plot",
                             choices = c("Dots", "Bars"))
               )),
        column(9,
               plotOutput('viewclonalbias', height = 700)
        )
      )
    })

    #======================================================================================================
    #SINGLECLONETRACKER TAB

    output$SingleTracker <- renderUI({
      if(is.null(thresholded_data()))
        return()

      singletrackerInput <- function(){
        print(barcodetrackR::single_clone_tracker(singletracker_data(),
                                                  n_clones = input$singletracker_Clones,
                                                  y_max = input$singletracker_yMax,
                                                  top_value = input$singletracker_Topvalue,
                                                  line_size = input$singletracker_Linesize,
                                                  text_size = input$singletracker_Textsize))
      }

      output$viewsingletracker<- renderPlot({
        singletrackerInput()
        height = 700
      })

      singletracker_data <- reactive({
        st_df <- thresholded_data()
        st_df <- st_df[st_df$GIVENNAME %in% input$singletracker_Samples,] #subset samples
        st_df$GIVENNAME <- factor(st_df$GIVENNAME, levels = input$singletracker_Samples)
        st_df <- st_df[order(st_df$GIVENNAME),]
        newcolnames <- st_df$GIVENNAME
        st_df$GIVENNAME <- NULL
        st_df$EXPERIMENT <- NULL
        st_df$CELLTYPE <- NULL
        st_df$MONTH <- NULL
        st_df$LOCATION <- NULL
        st_df$MISC <- NULL
        st_df <- data.frame(t(st_df))
        colnames(st_df) <- newcolnames
        return(st_df)
      })
      fluidRow(
        column(3,
               wellPanel(
                 selectInput("singletracker_Samples", label = "1. Which Samples to Use",
                             choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = TRUE),
                 numericInput("singletracker_Clones", "2. Enter Number of Top Clones: ", value = 100),
                 numericInput("singletracker_yMax", "3. Enter Y axis Limit: ", value = 1),
                 numericInput("singletracker_Topvalue", "4. Enter Clone Peak Choice: ", value = .5),
                 numericInput("singletracker_Linesize", "5. Enter Line Size: ", value = 2),
                 numericInput("singletracker_Textsize", "6. Enter Text Size: ", value = 10)
               )),
        column(9,
               plotOutput('viewsingletracker', height = 700)
        )
      )
    })




  }
)





