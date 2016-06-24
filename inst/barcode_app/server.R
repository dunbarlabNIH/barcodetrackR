options(shiny.maxRequestSize=1000*1024^2)

shinyServer(

  function(input, output) {

    #=========================================================================================================

    output$thresholdPanel <- renderUI({
      if (is.null(input$file1) | is.null(input$file2))
        return()

      column(3,
             wellPanel(
               strong("4. Press Button to Apply Threshold"),
               actionButton("threshybutton", "Apply Threshold")
             ))

    })

    thresholded_data <- eventReactive(input$threshybutton, {
      withProgress(message = "Applying Threshold", value = 0, {
        df <- barcodetrackR::threshold(read.delim(input$file1$datapath, row.names = 1), input$thresholdvalue)
        incProgress(0.5, detail = "Applying Threshold")
        tdf <- t(df)
        keyfile <- read.delim(input$file2$datapath, row.names = 1)
        tdf <- tdf[rownames(tdf) %in% rownames(keyfile),]
        keyfile = keyfile[match(rownames(tdf), rownames(keyfile)),]
        df <- cbind(keyfile, tdf)
      })
      df <- data.frame(df)
      return(df)
    }
    )


    current_threshold <- eventReactive(input$threshybutton,{
      return(paste("Current threshold applied is: ", input$thresholdvalue, sep = ""))
    })


    output$thresholdInfo <- renderUI({
      if (is.null(thresholded_data()))
        return()

      column(3, align = "center",
             wellPanel(
               strong(current_threshold())
             ))


    })
    #======================================================================================================

    #TABPANEL

    output$Panel <- renderUI({
      if (is.null(thresholded_data()))
        return()

      tabsetPanel(
        tabPanel("BCheatmap",
                 uiOutput("BCheatmap")),
        tabPanel("CorPlot",
                 uiOutput("CorPlot")),
        tabPanel("RadarChart",
                 uiOutput("RadarChart")),
        tabPanel("TopClonesBarChart",
                 uiOutput("TopClonesBarChart")),
        tabPanel("TopClonesTracker",
                 uiOutput("TopClonesTracker")),
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

    #HEATMAP TAB

    output$BCheatmap <- renderUI({

      if (is.null(thresholded_data()))
        return()


      #======================================================================================================

      BCheatmapInput <- function(){
        if(input$BCheatmap_cluster_tracker){
          print(BC_clusters(your_data = BCheatmap_data(),
                            names = paste0(colnames(BCheatmap_data()), "  "),
                            n_clones = input$BCheatmap_top_clones,
                            months = BCheatmap_months(),
                            your_title = input$BCheatmap_title,
                            log_transform = input$BCheatmap_log_transform,
                            log_choice = switch(as.character(input$BCheatmap_scale), "2" = 2, "e" = exp(1), "10" = 10, "100" = 100),
                            distance_method = input$BCheatmap_distance,
                            minkowski_power = input$BCmink_distance,
                            hclust_linkage = input$BCheatmap_hclust_linkage,
                            clusters = input$BCheatmap_clusters,
                            variable_log_min = input$BCheatmap_minimum,
                            percent_scale = input$BCheatmap_clusterscale
          ))
        } else {
          barcodetrackR::BCheatmap(your_data = BCheatmap_data(),
                                   names = paste0(colnames(BCheatmap_data()), "  "),
                                   n_clones = input$BCheatmap_top_clones,
                                   your_title = input$BCheatmap_title,
                                   grid = input$BCheatmap_grid, columnLabels = input$BCheatmap_labels,
                                   star_size = input$BCheatmap_starsize, log_transform = input$BCheatmap_log_transform,
                                   log_choice = switch(as.character(input$BCheatmap_scale), "2" = 2, "e" = exp(1), "10" = 10, "100" = 100),
                                   distance_method = input$BCheatmap_distance,
                                   minkowski_power = input$BCmink_distance,
                                   cellnote_option = input$BCheatmap_cellnote_option,
                                   hclust_linkage = input$BCheatmap_hclust_linkage,
                                   row_order = input$BCheatmap_row_order,
                                   clusters = input$BCheatmap_clusters,
                                   dendro = input$BCheatmap_dendrogram,
                                   variable_log_min = input$BCheatmap_minimum
          )
        }
      }

      output$downloadBCheatmapkey <- downloadHandler(
        filename = function() {paste(input$file1, "_BCheatmapkey.txt", sep = "")},
        content = function(file){
          write.table(barcodetrackR::BCheatmap(your_data = BCheatmap_data(),
                                               names = colnames(BCheatmap_data()),
                                               n_clones = input$BCheatmap_top_clones,
                                               log_transform = input$BCheatmap_log_transform,
                                               printtable = TRUE,
                                               table_option = input$BCheatmap_table_option,
                                               log_choice = switch(as.character(input$BCheatmap_scale), "2" = 2, "e" = exp(1), "10" = 10),
                                               distance_method = input$BCheatmap_distance,
                                               minkowski_power = input$BCmink_distance,
                                               hclust_linkage = input$BCheatmap_hclust_linkage,
                                               variable_log_min = input$BCheatmap_minimum
          ), file, sep = '\t', quote = FALSE)
        }
      )

      output$viewBCheatmap <- renderPlot({
        BCheatmapInput()
        height = 700
      })


      BCheatmap_data <- reactive({
        df <- thresholded_data()
        df <- df[df$GIVENNAME %in% input$BCheatmap_samples,] #subset samples
        df$GIVENNAME <- factor(df$GIVENNAME, levels = input$BCheatmap_samples)
        df <- df[order(df$GIVENNAME),]
        newcolnames <- df$GIVENNAME
        df$GIVENNAME <- NULL
        df$EXPERIMENT <- NULL
        df$CELLTYPE <- NULL
        df$MONTH <- NULL
        df$LOCATION <- NULL
        df$MISC <- NULL
        df <- data.frame(t(df))
        colnames(df) <- newcolnames
        return(df)

      })

      BCheatmap_months <- reactive({
        return(as.numeric(unlist(strsplit(input$BCheatmap_months, split = ','))))
      })

      fluidRow(
        column(3,
               wellPanel(
                 selectizeInput("BCheatmap_samples", label = "1. Which Samples to Use (in order)",
                                choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = TRUE),
                 textInput("BCheatmap_months", '(1.5). Enter months (in order) seperated by a comma: ', value = ""),
                 numericInput("BCheatmap_top_clones", "2. Number of top clones", value = 10),
                 textInput("BCheatmap_title", "3. Title for BCheatmap", value = ""),
                 strong("4. Options"),
                 checkboxInput("BCheatmap_grid", label = "Grid", value = TRUE),
                 checkboxInput("BCheatmap_log_transform", label = "Log-Transform", value = TRUE),
                 checkboxInput("BCheatmap_dendrogram", label = "Display Dendrogram", value = FALSE),
                 checkboxInput("BCheatmap_minimum", label = "Variable Minimum for Log", value = TRUE),
                 checkboxInput("BCheatmap_cluster_tracker", label = "Show cluster-tracker (CT)", value = FALSE),
                 checkboxInput("BCheatmap_clusterscale", label = "CT: percent scale", value = TRUE),
                 checkboxInput("BCheatmap_cluster_error", label = "CT: use std_error", value = TRUE),
                 selectInput("BCheatmap_scale", "4. Select Log",
                             choices = c("2", "e", "10", "100"),
                             selected = "e"),
                 selectInput("BCheatmap_distance", "5. Select Distance Metric/Function",
                             choices = sort(as.vector(unlist(summary(proxy::pr_DB)[1]))),
                             selected = "Euclidean"),
                 numericInput("BCmink_distance", "6. If Minkowski, choose Minkowski Power", value = 2, step = 1),
                 selectInput("BCheatmap_cellnote_option", "7. Select Cell Display Option",
                             choices = c("reads", "percents", "logs", "stars", "ranks"),
                             selected = "stars"),
                 selectInput("BCheatmap_hclust_linkage", "8. Select Clustering Linkage",
                             choices = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"),
                             selected = "complete"),
                 numericInput("BCheatmap_clusters", "9. Select number of clusters to cut", value = 0, step = 1, min = 1, max = 9),
                 selectInput("BCheatmap_row_order", "10. How to order rows", choices = c("hierarchical", "emergence"), selected = "hierarchical"),
                 numericInput("BCheatmap_labels", "10. Set Column Label Size", value = 3),
                 numericInput("BCheatmap_starsize", "11. Set Cell Label Size", value = 1.5),
                 selectInput("BCheatmap_table_option", "12. Select Format for Key (below)",
                             choices = c("logs", "reads", "percents", "ranks"),
                             selected = "percents"),
                 strong("13. Press button to download BCheatmap Key."),
                 br(),
                 downloadButton('downloadBCheatmapkey', 'BCheatmap_key')




               )
        ),
        column(9,
               plotOutput('viewBCheatmap', height = 900)
        )

      )



    })

    #======================================================================================================

    #CORPLOT TAB


    output$CorPlot <- renderUI({

      if (is.null(thresholded_data()))
        return()

      corplotInput <- function(){
        barcodetrackR::cor_plot(your_data = corplot_data(), names = colnames(corplot_data()),
                                thresh = input$corplot_thresh, your_title = input$corplot_Title,
                                method_corr = input$corplot_Method, labelsizes = input$corplot_Labels,
                                plottype = input$corplot_Type,
                                no_negatives = input$corplot_excludeneg,
                                show_grid = input$corplot_Grid, colorscale = input$corplot_Colors)
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
        height = 900
      })

      corplot_data <- reactive({
        cf <- thresholded_data()
        cf <- cf[cf$GIVENNAME %in% input$corplot_Samples,] #subset samples
        cf$GIVENNAME <- factor(cf$GIVENNAME, levels = input$corplot_Samples)
        cf <- cf[order(cf$GIVENNAME),]
        newcolnames <- cf$GIVENNAME
        cf$GIVENNAME <- NULL
        cf$EXPERIMENT <- NULL
        cf$CELLTYPE <- NULL
        cf$MONTH <- NULL
        cf$LOCATION <- NULL
        cf$MISC <- NULL
        cf <- data.frame(t(cf))
        colnames(cf) <- newcolnames
        return(cf)

      })

      #======================================================================================================


      fluidRow(
        column(3,
               wellPanel(
                 selectizeInput("corplot_Samples", label = "1. Which Samples to Use (order matters)",
                                choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = TRUE),
                 numericInput("corplot_thresh", "2. Select Threshold (for Corrplot)", value = 0),
                 textInput("corplot_Title", "3. Title for Corplot", value = ""),
                 selectInput("corplot_Type", '4. Choose Type of Plot', choices = c("circle", "square", "ellipse", "color", "number", "shade", "pie"), selected = "square"),
                 strong("5. Exclude Negatives?"),
                 checkboxInput("corplot_excludeneg", "", value = FALSE),
                 strong("6. Grid ON/OFF"),
                 checkboxInput("corplot_Grid", "", value = TRUE),
                 selectInput("corplot_Method", "7. Chooose Correlation Method", choices = c("pearson", "kendall", "spearman", "manhattan"), selected = "pearson"),
                 selectInput("corplot_Colors", "8. Choose Color Scale", choices = c("default", "rainbow", "white_heat"), selected = "default"),
                 numericInput("corplot_Labels", "9. Set Label Size", value = 2),
                 strong("10. Press button to downlaod CorPlot Files as .zip"),
                 br(),
                 downloadButton('downloadcorplotzip', 'corr_plot.zip')
               )
        ),

        column(9,
               plotOutput('viewcorplot', height = 1000)
        )

      )





    })


    #======================================================================================================




    output$RadarChart <- renderUI({


      radarchart_data <- reactive({

        rf <- thresholded_data()
        rf <- rf[rf$GIVENNAME %in% input$radar_Samples,] #subset samples
        rf$GIVENNAME <- factor(rf$GIVENNAME, levels = input$radar_Samples)
        rf <- rf[order(rf$GIVENNAME),]
        newcolnames <- rf$GIVENNAME
        rf$GIVENNAME <- NULL
        rf$EXPERIMENT <- NULL
        rf$CELLTYPE <- NULL
        rf$MONTH <- NULL
        rf$LOCATION <- NULL
        rf$MISC <- NULL
        rf <- data.frame(t(rf))
        colnames(rf) <- newcolnames
        return(rf)

      })

      radarInput <- function(){
        barcodetrackR::radartopclones(your_data = radarchart_data(),
                                      columnChoice_Name = input$radar_Choice,
                                      n_clones = input$radar_Clones, your_title = input$radar_Title,
                                      labelsize = input$radar_labelsize)
      }
      output$viewradarchart <- renderPlot({
        radarInput()
      })


      fluidRow(
        column(3,
               wellPanel(
                 selectInput("radar_Samples", label = "1. Which Samples to Use",
                             choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = TRUE),

                 selectizeInput("radar_Choice", label = "2. Which Sample for Top Clones: ",
                                choices = as.vector(unique(thresholded_data()$GIVENNAME))),

                 numericInput("radar_Clones", "3. Select Number of Clones", value = 100, step = 10),

                 textInput("radar_Title", "4. Title", value = ""),

                 numericInput("radar_labelsize", "5. Enter Label Size: ", value = 1)



               )),




        column(9,
               plotOutput('viewradarchart', height = 700)
        )

      )
    })





    #======================================================================================================

    #TOPCLONESBARCHART TAB

    output$TopClonesBarChart <- renderUI({

      if (is.null(thresholded_data()))
        return()

      topclonesbarchartInput <- function(){
        barcodetrackR::topclones_barchart(your_data = topclonesbarchart_data(),
                                          top_clones_choice = topclonesbarchart_Choice(),
                                          n_clones = input$topclonesbarchart_Clones,
                                          text_size = input$topclonesbarchart_Textsize,
                                          y_limit = input$topclonesbarchart_yLim,
                                          other_color = input$topclonesbarchart_Othercolor,
                                          your_title = input$topclonesbarchart_Title)
      }


      output$viewtopclonesbarchart <- renderPlot({
        print(topclonesbarchartInput())
        height = 700
      })

      topclonesbarchart_data <- reactive({
        ef <- thresholded_data()
        ef <- ef[ef$GIVENNAME %in% input$topclonesbarchart_Samples,] #subset samples
        ef$GIVENNAME <- factor(ef$GIVENNAME, levels = input$topclonesbarchart_Samples)
        ef <- ef[order(ef$GIVENNAME),]
        newcolnames <- ef$GIVENNAME
        ef$GIVENNAME <- NULL
        ef$EXPERIMENT <- NULL
        ef$CELLTYPE <- NULL
        ef$MONTH <- NULL
        ef$LOCATION <- NULL
        ef$MISC <- NULL
        ef <- data.frame(t(ef))
        colnames(ef) <- newcolnames
        return(ef)
      })

      topclonesbarchart_Choice <- reactive({
        tcbchoice <- thresholded_data()
        tcbchoice <- tcbchoice[tcbchoice$GIVENNAME == input$topclonesbarchart_Selected,,drop = FALSE]
        newcolnames <- tcbchoice$GIVENNAME
        tcbchoice$GIVENNAME <- NULL
        tcbchoice$EXPERIMENT <- NULL
        tcbchoice$CELLTYPE <- NULL
        tcbchoice$MONTH <- NULL
        tcbchoice$LOCATION <- NULL
        tcbchoice$MISC <- NULL
        tcbchoice <- data.frame(t(tcbchoice))
        colnames(tcbchoice) <- newcolnames
        return(tcbchoice)
      })


      fluidRow(
        column(3,
               wellPanel(
                 selectInput("topclonesbarchart_Samples", label = "1. Which Samples to Use",
                             choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = TRUE),
                 selectInput("topclonesbarchart_Selected", label = "2. Which Sample to Use for Top Clones",
                             choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = FALSE),
                 numericInput("topclonesbarchart_Clones", "3. Select Number of Clones", value = 50),
                 numericInput("topclonesbarchart_Textsize", "4. Enter Text Size: ", value = 15),
                 numericInput("topclonesbarchart_yLim", "5. Enter Y Limit: ", value = 100),
                 selectInput("topclonesbarchart_Othercolor", "6. Choose Other Color",
                             choices = c("grey", "black", "white"), selected = "black"),
                 textInput("topclonesbarchart_Title", "7. Title", value = "Bar Chart")



               )),



        column(9,
               plotOutput('viewtopclonesbarchart', height = 700)
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



        column(9,
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






  }
)





