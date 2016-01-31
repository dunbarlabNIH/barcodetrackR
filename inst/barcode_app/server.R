options(shiny.maxRequestSize=1000*1024^2)

shinyServer(
  
  function(input, output) {
    
    
    #=================================================================================``insta========================
    
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
        tabPanel("TopClonesContribt",
                 uiOutput("TopClonesContrib")),
        tabPanel("HematoPer",
                 uiOutput("HematoPer")),
        tabPanel("Diversity & Richness",
                 uiOutput("Diversity")),
        tabPanel("Scatter Plot",
                 uiOutput("ScatterPlot")),
        tabPanel("Tree Map",
                 uiOutput("TreeMap")),
        tabPanel("BBHM",
                 uiOutput("BBHM"))
      )
      
    })
    
    #======================================================================================================
    
    #HEATMAP TAB
    
    output$BCheatmap <- renderUI({
      
      if (is.null(thresholded_data()))
        return()
      
      
      #======================================================================================================
      
      BCheatmapInput <- function(){
        barcodetrackR::BCheatmap(your_data = BCheatmap_data(),
                                 names = colnames(BCheatmap_data()),
                                 n_clones = input$BCheatmap_top_clones,
                                 your_title = input$BCheatmap_title,
                                 grid = input$BCheatmap_grid, columnLabels = input$BCheatmap_labels,
                                 star_size = input$BCheatmap_starsize, log_transform = input$BCheatmap_log_transform,
                                 log_choice = switch(as.character(input$BCheatmap_scale), "2" = 2, "e" = exp(1), "10" = 10, "100" = 100),
                                 distance_method = input$BCheatmap_distance,
                                 minkowski_power = input$BCmink_distance,
                                 cellnote_option = input$BCheatmap_cellnote_option,
                                 hclust_linkage = input$BCheatmap_hclust_linkage
        )
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
                                               hclust_linkage = input$BCheatmap_hclust_linkage
          ), file, sep = '\t', quote = FALSE)
        }
      )
      
      output$viewBCheatmap <- renderPlot({
        BCheatmapInput()
        height = 900
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
      
      fluidRow(
        column(3,
               wellPanel(
                 selectizeInput("BCheatmap_samples", label = "1. Which Samples to Use (in order)",
                                choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = TRUE),
                 numericInput("BCheatmap_top_clones", "2. Number of top clones", value = 10),
                 textInput("BCheatmap_title", "2. Title for BCheatmap", value = ""),
                 strong("3. Grid ON/OFF"),
                 checkboxInput("BCheatmap_grid", "", value = TRUE),
                 strong("4. Log Transform?"),
                 checkboxInput("BCheatmap_log_transform", "", value = TRUE),
                 selectInput("BCheatmap_scale", "4. Select Log",
                             choices = c("2", "e", "10", "100"),
                             selected = "e"),
                 selectInput("BCheatmap_distance", "5. Select Distance Metric",
                             choices = c("euclidean", "maximum", "manhattan", "canberra", "binary","minkowski"),
                             selected = "euclidean"),
                 numericInput("BCmink_distance", "6. If Minkowski, choose Minkowski Power", value = 2, step = 1),
                 selectInput("BCheatmap_cellnote_option", "7. Select Cell Display Option",
                             choices = c("reads", "percents", "logs", "stars"),
                             selected = "stars"),
                 selectInput("BCheatmap_hclust_linkage", "8. Select Clustering Linkage",
                             choices = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"),
                             selected = "complete"),
                 numericInput("BCheatmap_labels", "9. Set Column Label Size", value = 1.5),
                 numericInput("BCheatmap_starsize", "10. Set Cell Label Size", value = 1.5),
                 selectInput("BCheatmap_table_option", "11. Select Format for Key (below)",
                             choices = c("logs", "reads", "percents"),
                             selected = "percents"),
                 strong("12. Press button to download BCheatmap Key."),
                 br(),
                 downloadButton('downloadBCheatmapkey', 'BCheatmap_key')
                 
                 
               )
        ),
        column(9,
               plotOutput('viewBCheatmap', height = 1000)
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
                 selectInput("corplot_Method", "7. Chooose Correlation Method", choices = c("pearson", "kendall", "spearman"), selected = "pearson"),
                 selectInput("corplot_Colors", "8. Choose Color Scale", choices = c("default", "rainbow", "white_heat"), selected = "default"),
                 numericInput("corplot_Labels", "9. Set Label Size", value = 0.5),
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
    
    #TOPCLONESCONTRIB TAB
    
    output$TopClonesContrib <- renderUI({
      
      if (is.null(thresholded_data()))
        return()
      
      
      topclonescontribInput <- function(){
        barcodetrackR::topclonescontrib(your_data = topclonescontrib_data(),
                                        n_clones = input$topclonescontrib_Clones, linesize = input$topclonescontrib_Linesize,
                                        pointsize = input$topclonescontrib_Pointsize, your_title = input$topclonescontrib_Title)
      }
      
      
      output$viewtopclones <- renderPlot({
        print(topclonescontribInput())
        height = 700
      })
      
      topclonescontrib_data <- reactive({
        
        ef <- thresholded_data()
        
        ef <- ef[ef$GIVENNAME %in% input$topclonescontrib_Samples,] #subset samples
        ef$GIVENNAME <- factor(ef$GIVENNAME, levels = input$topclonescontrib_Samples)
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
      
      
      fluidRow(
        column(3,
               wellPanel(
                 selectInput("topclonescontrib_Samples", label = "1. Which Samples to Use",
                             choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = TRUE),
                 
                 numericInput("topclonescontrib_Clones", "2. Select Number of Clones", value = 50),
                 
                 numericInput("topclonescontrib_Linesize", "3. Select Line Size", value = 2),
                 
                 numericInput("topclonescontrib_Pointsize", "4. Select Point Size", value = 3),
                 
                 textInput("topclonescontrib_Title", "5. Title for Plot", value = "")
                 
                 
                 
               )),  
        
        
        
        column(9,
               plotOutput('viewtopclones', height = 700)
        )
        
      )
      
      
      
    })
    
    
    #======================================================================================================
    #HEMATOPER TAB
    
    output$HematoPer <- renderUI({
      if(is.null(thresholded_data()))
        return()
      
      hematoperInput <- function(){
        print(barcodetrackR::hematoper(hematoper_data(), months = hematoper_Months(), n_clones = input$hematoper_Clones, scale_percent = input$hematoper_Scale, linesize = input$hematoper_Linesize))
        
      }
      
      output$viewhematoper <- renderPlot({
        hematoperInput()
        height = 700
      })
      
      hematoper_data <- reactive({
        ff <- thresholded_data()
        
        ff <- ff[ff$GIVENNAME %in% input$hematoper_Samples,] #subset samples
        ff$GIVENNAME <- factor(ff$GIVENNAME, levels = input$hematoper_Samples)
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
      
      hematoper_Months <- reactive({
        return(as.numeric(unlist(strsplit(input$hematoper_Months, split = ','))))
      })
      
      fluidRow(
        column(3,
               wellPanel(
                 selectInput("hematoper_Samples", label = "1. Which Samples to Use",
                             choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = TRUE),
                 
                 textInput("hematoper_Months", '2. Enter months (in order) seperated by a comma: ', value = ""),
                 
                 numericInput("hematoper_Clones", "3. Select Number of Clones", value = 10),
                 
                 numericInput("hematoper_Linesize", "4. Enter Line Size: ", value = 1),
                 
                 strong("5. Scale as Percent?"),
                 
                 checkboxInput("hematoper_Scale", "", value = TRUE)
                 
               )),  
        
        
        
        column(9,
               plotOutput('viewhematoper', height = 700)
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
    
    
    
  }
)



