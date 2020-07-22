options(shiny.maxRequestSize=1000*1024^2)
library(SummarizedExperiment)

shinyServer(

  function(input, output, session) {

    #=========================================================================================================

    output$thresholdPanel <- renderUI({
      if (is.null(input$file1) | is.null(input$file2)) #| is.null(input$file3))
        return()
      actionButton("threshybutton", "Load Files and Apply Threshold", width="100%")


    })

    thresholded_data <- eventReactive(input$threshybutton, {
      withProgress(message = "Loading files and applying threshold", value = 0, {
        # your_data <- barcodetrackR::threshold(read.delim(input$file1$datapath, row.names = 1), input$thresholdvalue)
        your_data <- read.delim(input$file1$datapath,  row.names =1)
        metadata <- read.delim(input$file2$datapath)

        # incProgress(0.5, detail = "Loading metadata/readme and applying Threshold")
        # t_your_data <- t(your_data)
        # metadata <- read.delim(input$file2$datapath, stringsAsFactors = FALSE)
        if(!(all(c("SAMPLENAME") %in% colnames(metadata)))){
          stop("metadata missing SAMPLENAME column")
        }
        # metadata <- metadata[,c("SAMPLENAME")]
        if(any(duplicated(metadata$SAMPLENAME))){
          stop("metadata contains duplicate SAMPLENAME")
        }
        # if(any(duplicated(metadata$GIVENNAME))){
        #   stop("metadata contains duplicate GIVENNAME")
        # }
        if(!(all(colnames(your_data) %in% metadata$SAMPLENAME))){
          stop("Column in data is not a SAMPLENAME in metadata")
        }
        if(!(all(metadata$SAMPLENAME %in% colnames(your_data)))){
          stop("SAMPLENAME in metadata is not a column in data")
        }
        if(length(setdiff(metadata$SAMPLENAME, colnames(your_data))) != 0){
          print(setdiff(metadata$SAMPLENAME, colnames(your_data)))
          stop("Number of samples in metadata differs from number of samples in data")
        }
        metadata <- metadata[match(colnames(your_data), metadata$SAMPLENAME),]
        your_SE <- barcodetrackR::create_SE(your_data = your_data,
                                            meta_data  = metadata,
                                            threshold = input$thresholdvalue)
        # t_your_data <- data.frame(GIVENNAME = metadata$GIVENNAME, t_your_data)
      })
      return(your_SE)
    }
    )


    current_threshold <- eventReactive(input$threshybutton,{
      return(paste("Current threshold applied is:\n", paste0(input$thresholdvalue*100, "%"), sep = ""))
    })


    output$thresholdInfo <- renderUI({
      if (is.null(thresholded_data()))
        return()
      strong(current_threshold())
    })

    # readme_data <- eventReactive(input$threshybutton,{
    #   if(is.null(thresholded_data()) | is.null(input$file3$datapath)){
    #     return()
    #   }
    #   readme <- read.delim(input$file3$datapath, stringsAsFactors = FALSE, comment.char = "#")
    #   if(!all(c("FILENAME", "MAPPED", "READS") %in% colnames(readme))){
    #     stop("Uploaded readme must contain FILENAME, MAPPED, and READS columns")
    #   }
    #   readme <- readme[, c("FILENAME", "MAPPED", "READS")]
    #   colnames(readme) <- c("FILENAME", "READS_WITH_LIBID", "RAW_READS")
    #   metadata <- read.delim(input$file2$datapath, stringsAsFactors = FALSE)
    #   readme$GIVENNAME <- metadata$GIVENNAME[match(readme$FILENAME, metadata$FILENAME)]
    #   readme$THRESHOLDED_READS <- rowSums(thresholded_data()[match(rownames(thresholded_data()), readme$FILENAME),-which(colnames(thresholded_data()) == "GIVENNAME")])
    #   readme <- readme[,c("FILENAME", "GIVENNAME", "READS_WITH_LIBID", "THRESHOLDED_READS", "RAW_READS")]
    #   readme$READS_WITH_LIBID_PERCENT <- round(readme$READS_WITH_LIBID/readme$RAW_READS * 100, digits = 2)
    #   readme$THRESHOLDED_READS_PERCENT <- round(readme$THRESHOLDED_READS/readme$RAW_READS * 100, digits = 2)
    #   return(readme)
    # })



    #======================================================================================================

    #TABPANEL

    output$Panel <- renderUI({

      if (is.null(thresholded_data()))
        return()

      tabsetPanel(
        tabPanel("Descriptive Statistics",
                  uiOutput("DataStatistics")),
        tabPanel("Heatmap",
                 uiOutput("Heatmap")),
        tabPanel("Correlation Plot",
                 uiOutput("CorPlot")),
        tabPanel("Clonal Diversity",
                 uiOutput("ClonalDiversity")),
        tabPanel("Clonal Contribution",
                 uiOutput("ClonalContribution")),
        tabPanel("Chord Diagram",
                 uiOutput("ChordDiagram")),
        tabPanel("Ridge Plot",
                 uiOutput("RidgePlot"))

        # tabPanel("Scatter Plot",
        #          uiOutput("ScatterPlot")),
        # tabPanel("Binary Heatmap",
        #          uiOutput("BinaryHeatmap")),
        # tabPanel("TernPlot",
        #          uiOutput("TernPlot")),
        # tabPanel("RankAbundance",
        #          uiOutput("RankAbundance")),
        # tabPanel("ClonalBias",
        #          uiOutput("ClonalBias"))
      )
    })


    #======================================================================================================

    #DATASTATISTICS TAB
    output$DataStatistics <- renderUI({
      
    if (is.null(thresholded_data()))
      return()
    
      stat_histInput <- function(){
        print(barcodetrackR::stat_hist(your_SE = thresholded_data(),
                                       data_choice = input$stat_hist_data_choice, 
                                       assay_choice = input$stat_hist_assay_choice, 
                                       sample_select = input$stat_hist_sample_select,
                                       metadata_stat = input$stat_hist_metadata_stat,
                                       group_by = input$stat_hist_group_by,
                                       scale_all_y = input$stat_hist_scale_all_y,
                                       y_log_axis = input$stat_hist_y_log_axis,
                                       text_size = input$stat_hist_text_size,
                                       n_bins = input$stat_hist_n_bins,
                                       n_cols = input$stat_hist_n_cols,
                                       alpha = input$stat_hist_alpha
        ))}
      
      output$viewDataStatistics <- renderPlot({
        stat_histInput()
      })
      
      fluidRow(
        column(3,
               wellPanel(
                 # div(style="display:inline-block; height:85px;",fileInput("Heatmap_uploaded_samples", "1. Upload Prepared Sample List or Input Samples")),
                 selectInput("stat_hist_data_choice", "Select Histogram Display Option",
                             choices = c("barcode_stats", "aggregate_stats"),
                             selected = "barcode_stats"),
                 selectizeInput("stat_hist_sample_select", label = "Samples", choices = as.vector(unique(thresholded_data()$SAMPLENAME)), multiple = TRUE),
                 selectInput("stat_hist_assay_choice", label = "Choose Assay", choices = names(SummarizedExperiment::assays(thresholded_data())), multiple = FALSE),
                 numericInput("stat_hist_n_cols", "Number of columns", value = 2),
                 # br(),
                 selectizeInput("stat_hist_metadata_stat", label = "Choose Metadata", choices = colnames(SummarizedExperiment::colData(thresholded_data()))[unlist(lapply(SummarizedExperiment::colData(thresholded_data()), is.numeric))], multiple = FALSE),
                 selectInput("stat_hist_group_by", label = "Group by", choices = colnames(SummarizedExperiment::colData(thresholded_data())), multiple = FALSE, selected = FALSE, selectize = FALSE, size = 6), # somehow adding selectize = FALSE and size = some number allows for null default
                 # br(),
                 strong("Options"),
                 checkboxInput("stat_hist_scale_all_y", label = "Scale all y", value = FALSE),
                 checkboxInput("stat_hist_y_log_axis", label = "Y log axis", value = FALSE),
                 numericInput("stat_hist_text_size", "Text size", value = 12),
                 numericInput("stat_hist_n_bins", "Number of bins", value = 10),
                 numericInput("stat_hist_alpha", "Alpha", value = 0.5)
               )
        ),
        column(8,
               plotOutput('viewDataStatistics', height = 600)
        )
      )
    })                                
                                       
              
    #======================================================================================================
    #HEATMAP TAB

    output$Heatmap <- renderUI({

      if (is.null(thresholded_data()))
        return()

      HeatmapInput <- function(){
        print(barcodetrackR::barcode_ggheatmap(your_SE = Heatmap_data(),
                                               plot_labels = NULL,
                                               n_clones = input$Heatmap_top_clones,
                                               cellnote_assay = input$Heatmap_cellnote_assay,
                                               your_title = paste0("\n",input$Heatmap_title),
                                               grid = input$Heatmap_grid,
                                               label_size = input$Heatmap_labels,
                                               dendro = input$Heatmap_dendrogram,
                                               cellnote_size = input$Heatmap_starsize,
                                               distance_method = input$Heatmap_distance,
                                               minkowski_power = input$Heatmap_mink_distance,
                                               hclust_linkage = input$Heatmap_hclust_linkage,
                                               row_order = input$Heatmap_row_order,
                                               clusters = input$Heatmap_clusters,
                                               percent_scale = c(0, 0.000025, 0.001, 0.01, 0.1, 1),
                                               color_scale = c("#4575B4", "#4575B4", "lightblue", "#fefeb9", "#D73027", "red4")
                                               # log_transform = input$Heatmap_log_transform,
                                               # log_choice = switch(as.character(input$Heatmap_scale), "2" = 2, "e" = exp(1), "10" = 10, "100" = 100),
                                               ))}

      # HeatmapREADMEtable <- function(){
      #   your_table <- readme_data()
      #   your_table <- your_table[match(input$Heatmap_samples, your_table$GIVENNAME),]
      #   temp_names <- your_table$GIVENNAME
      #   temp_clones <- colSums(Heatmap_data() > 0)
      #   your_table <- data.frame(THRESHOLDED_READS = your_table$THRESHOLDED_READS, TOTAL_READS = your_table$RAW_READS, THRESHOLDED_READS_PERCENT = your_table$THRESHOLDED_READS_PERCENT, N_CLONES = temp_clones)
      #   rownames(your_table) <- temp_names
      #   return(cbind(rownames(t(your_table)),t(your_table)))
      # }


      # output$heatmap_datatable <- renderDataTable({HeatmapREADMEtable()},
      #                                             options = list(paging = FALSE, searching = FALSE))
      # 
      # 
      # output$downloadHeatmapkey <- downloadHandler(
      #   filename = function() {paste(input$file1, "_heatmapkey.txt", sep = "")},
      #   content = function(file){
      #     write.table(barcodetrackR::barcode_ggheatmap(your_SE = Heatmap_data(),
      #                                                  plot_labels = NULL,
      #                                                  n_clones = input$Heatmap_top_clones,
      #                                                  log_transform = input$Heatmap_log_transform,
      #                                                  printtable = TRUE,
      #                                                  table_option = input$Heatmap_table_option,
      #                                                  log_choice = switch(as.character(input$Heatmap_scale), "2" = 2, "e" = exp(1), "10" = 10),
      #                                                  distance_method = input$Heatmap_distance,
      #                                                  minkowski_power = input$Heatmap_mink_distance,
      #                                                  hclust_linkage = input$Heatmap_hclust_linkage), file, sep = '\t', quote = FALSE)
      #   }
      # )
      # 
      # output$downloadHeatmapSTATS <- downloadHandler(
      #   filename = function() {paste(input$file1, "_heatmapSTATS.txt", sep = "")},
      #   content = function(file){
      #     write.table(HeatmapREADMEtable(), file, sep = '\t', quote = FALSE)
      #   }
      # )
      
      
      output$viewHeatmap <- renderPlot({
        HeatmapInput()
      })
 
      Heatmap_data <- reactive({
        # df <- thresholded_data()
        se <- thresholded_data()
        se <- se[,se$SAMPLENAME %in% input$Heatmap_samples] # subset samples
        se$SAMPLENAME <- factor(se$SAMPLENAME, levels = input$Heatmap_samples)
        se <- se[,order(se$SAMPLENAME)]
        # df <- df[df$GIVENNAME %in% input$Heatmap_samples,] #subset samples
        # df$GIVENNAME <- factor(df$GIVENNAME, levels = input$Heatmap_samples)
        #df <- df[order(df$GIVENNAME),]
        # newcolnames <- df$GIVENNAME
        # df$GIVENNAME <- NULL
        # df <- data.frame(t(df))
        # colnames(df) <- newcolnames
        return(se)

      })

      # output$Heatmap_download_samples <- downloadHandler(
      #   filename = function() {paste("_samplelist.txt", sep = "\t")},
      #   content = function(file){
      #     write.table(input$Heatmap_samples, file, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
      #   }
      # )

      Heatmap_uploaded_samples <- reactive({
        samples <- as.vector(t(read.delim(input$Heatmap_uploaded_samples$datapath, stringsAsFactors = FALSE, header = FALSE)))
        return(samples)
      })

      observeEvent(input$Heatmap_uploaded_samples, {updateSelectizeInput(session, inputId = 'Heatmap_samples', selected = Heatmap_uploaded_samples())})

      fluidRow(
        column(3,
               wellPanel(
                 div(style="display:inline-block; height:85px;",fileInput("Heatmap_uploaded_samples", "Upload Prepared Sample List or Input Samples")),
                 selectizeInput("Heatmap_samples", label = NULL, choices = as.vector(unique(thresholded_data()$SAMPLENAME)), multiple = TRUE),
                 # downloadButton("Heatmap_download_samples", "Download This Sample List"),
                 numericInput("Heatmap_top_clones", "Number of top clones", value = 10),
                 strong("Options"),
                 checkboxInput("Heatmap_grid", label = "Grid", value = TRUE),
                 # checkboxInput("Heatmap_log_transform", label = "Log-Transform", value = TRUE),
                 checkboxInput("Heatmap_dendrogram", label = "Display Dendrogram", value = FALSE),
                 selectInput("Heatmap_cellnote_assay", "Select Cell Display Option",
                             choices = c("stars", "reads", "percentages"),
                             selected = "stars"),
                 numericInput("Heatmap_labels", "Set Column Label Size", value = 14),
                 numericInput("Heatmap_starsize", "Set Cell Label Size", value = 6),
                 textInput("Heatmap_title", "Title for Heatmap", value = ""),
                 selectInput("Heatmap_distance", "Select Distance Metric/Function",
                             choices = sort(as.vector(unlist(summary(proxy::pr_DB)[1]))),
                             selected = "Euclidean"),
                 numericInput("Heatmap_mink_distance", "If Minkowski, choose Minkowski Power", value = 2, step = 1),
                 selectInput("Heatmap_hclust_linkage", "Select Clustering Linkage",
                             choices = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"),
                             selected = "complete"),
                 # selectInput("Heatmap_scale", "11. Select Log",
                 #             choices = c("2", "e", "10", "100"),
                 #             selected = "e"),
                 numericInput("Heatmap_clusters", "Select number of clusters to show", value = 0, step = 1, min = 1, max = 11),
                 selectInput("Heatmap_row_order", "How to order rows", choices = c("hierarchical", "emergence"), selected = "hierarchical"),
                 # selectInput("Heatmap_table_option", "14. Select Format for Key (below)",
                 #             choices = c("logs", "reads", "percents", "ranks"),
                 #             selected = "percents"),
                 # strong("15. Press button to download Heatmap Key."),
                 # br(),
                 # downloadButton('downloadHeatmapkey', 'Heatmap_key'),
                 # br(),
                 # strong("16. Press button to download Heatmap Stats."),
                 # br(),
                 # downloadButton('downloadHeatmapSTATS', 'Heatmap_STATS')
               )
        ),
        column(8,
               plotOutput('viewHeatmap', height = 800),
            #   dataTableOutput('heatmap_datatable')
        )
      )
    })

    
    #======================================================================================================
    #CORPLOT TAB

    output$CorPlot <- renderUI({


      if (is.null(thresholded_data()))
        return()

     

      corplotInput <- function(){
        barcodetrackR::cor_plot(your_SE = corplot_data(),
                                # thresh = input$corplot_thresh,
                                your_title = input$corplot_Title,
                                method_corr = input$corplot_Method,
                                label_size = input$corplot_Labels,
                                plot_type = input$corplot_Type,
                                no_negatives = input$corplot_excludeneg,
                                grid = input$corplot_Grid,
                               # color_scale = input$corplot_Colors,
                                number_size = input$corplot_number_size,
                                point_scale = input$corplot_point_scale)
      }

      output$downloadcorplotzip <- downloadHandler(
        filename = function() {"corplot_files.zip"},
        content = function(file){
          tmpdir <- tempdir()
          setwd(tmpdir)
          listoffiles <- barcodetrackR::cor_plot(your_data = corplot_data(),
                                                 # thresh = input$corplot_thresh, 
                                                 your_title = input$corplot_Title,
                                                 method_corr = input$corplot_Method,
                                                 label_size = input$corplot_Labels,
                                                 plot_type = input$corplot_Type,
                                                 return_table = TRUE)
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
        se <- thresholded_data()
        se <- se[,se$SAMPLENAME %in% input$corplot_Samples] # subset samples
        se$SAMPLENAME <- factor(se$SAMPLENAME, levels = input$corplot_Samples)
        se <- se[,order(se$SAMPLENAME)]
        return(se)
      })

      #======================================================================================================
      observeEvent(input$corplot_Copier, {updateSelectizeInput(session, inputId = 'corplot_Samples', selected = input$Heatmap_samples)})
      

      fluidRow(
        column(3,
               wellPanel(
                 selectizeInput("corplot_Samples", "Samples", choices = as.vector(unique(thresholded_data()$SAMPLENAME)), multiple = TRUE),
                 actionButton("corplot_Copier", label = "Copy Samples from heatmap"),
                 br(),
                 br(),
                 #numericInput("corplot_thresh", "2. Reads threshold", value = 0),
                 selectInput("corplot_Method", "Chooose Correlation Method", choices = c("pearson", "kendall", "spearman", "manhattan"), selected = "pearson"),
                 selectInput("corplot_Type", 'Choose Type of Plot', choices = c("color", "circle","number"), selected = "color"), 
                 textInput("corplot_Title", "Title", value = ""),
                 strong("Options"),
                 checkboxInput("corplot_excludeneg", "Exclude negatives", value = FALSE),
                 checkboxInput("corplot_Grid", "Grid", value = TRUE),
                # selectInput("corplot_Colors", "Choose Color Scale", choices = c("default", "rainbow", "white_heat"), selected = "default"),
                 numericInput("corplot_Labels", "Set Label Size", value = 12),
                 numericInput("corplot_number_size", "Set Number Size", value = 4),
                 numericInput("corplot_point_scale", "Set Point Scale", value = 4),
                 strong("Press button to downlaod CorPlot Files as .zip"),
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
    # Clonal Contribution TAB
    output$ClonalContribution <- renderUI({
      
      if (is.null(thresholded_data()))
        return()
      
      ClonalContributionInput <- function(){
        print(barcodetrackR::clonal_contribution(your_SE = thresholded_data(),
                                                 graph_type = input$cc_graph_type,
                                                 filter_by = input$cc_filter_by,
                                                 filter_selection = input$cc_filter_selection,
                                                 SAMPLENAME_choice = input$cc_samplename_choice,
                                                 plot_over = input$cc_plot_over,
                                                 # plot_over_display_choices = input$cc_plot_over_choices,
                                                 n_clones = input$cc_n_clones,
                                                 keep_numeric = input$cc_keep_numeric,
                                                 plot_non_selected = input$cc_plot_non_selected,
                                                 linesize = input$cc_linesize,
                                                 text_size = input$cc_text_size,
                                                 your_title = input$cc_your_title,
                                                 y_limit = input$cc_y_limit
        ))
      }
      
      
      output$view_cc_plot <- renderPlot({
        ClonalContributionInput()
      })
      
      
      # Helper function
      # filter_data_available <- reactive({
      #   filter_data_available <- unique(colData(thresholded_data())[,input$cc_filter_by])
      # })
      
      observeEvent(input$cc_filter_by, {updateSelectizeInput(session,
                                                             inputId = 'cc_filter_selection',
                                                             choices = unique(SummarizedExperiment::colData(thresholded_data())[,input$cc_filter_by]))})
      
      # observeEvent(input$cc_plot_over, {updateSelectizeInput(session,
      #                                                        inputId = 'cc_plot_over_choices',
      #                                                        choices = unique(colData(thresholded_data())[,input$cc_plot_over]))})

      fluidRow(
        column(3,
               wellPanel(
                 selectInput("cc_graph_type", "Chooose Graph Type", choices = c("bar", "line"), selected = "bar"),
                 selectInput("cc_filter_by", label = "Filter By", choices = colnames(SummarizedExperiment::colData(thresholded_data())), multiple = FALSE), 
                 selectInput("cc_filter_selection", label = "Filter Selection", choices = c(""), multiple = FALSE), 
                 selectizeInput("cc_samplename_choice", "Sample name to color by", choices = as.vector(unique(thresholded_data()$SAMPLENAME)), multiple = FALSE),
                 selectInput("cc_plot_over", label = "Plot Over", choices = colnames(SummarizedExperiment::colData(thresholded_data())), multiple = FALSE), 
                # selectizeInput("cc_plot_over_choices", label = "Choose values to plot over", choices = c(""), multiple = TRUE), 
                 numericInput("cc_n_clones", "Number of clones to color", value = 10),
                 textInput("cc_your_title", "Title", value = ""),
                 strong("Options"),
                 checkboxInput("cc_keep_numeric", "Keep x-axis numeric", value = TRUE),
                 checkboxInput("cc_plot_non_selected", "Plot non-selected", value = FALSE),
                 numericInput("cc_linesize", "Line size", value = 0.2),
                 numericInput("cc_text_size", "Text size", value = 16),
                 numericInput("cc_y_limit", "Manually set y limit", value = NULL)
               )
        ),
        
        column(8,
               plotOutput('view_cc_plot', height = 500)
        )
      )

    })
    

    #======================================================================================================
    # Clonal Diversity TAB
    output$ClonalDiversity <- renderUI({
      
      if (is.null(thresholded_data()))
        return()
      
      ClonalDiversityInput <- function(){
        print(barcodetrackR::clonal_diversity(your_SE = thresholded_data(),
                                              index_type = input$div_index,
                                              group_by = input$div_group_by,
                                              plot_over = input$div_plot_over,
                                              line_size = input$div_line_size,
                                              point_size = input$div_point_size,
                                              text_size = input$div_text_size,
                                              your_title = input$div_your_title
        ))
      }
      
      
      output$view_div_plot <- renderPlot({
        ClonalDiversityInput()
      })
      
      
      fluidRow(
        column(3,
               wellPanel(
                 selectInput("div_index", "Chooose Diversity Index",
                             choices = c("count","shannon","shannon_count","simpson", "invsimpson"), selected = "shannon"),
                 selectInput("div_group_by", label = "Group By", choices = colnames(SummarizedExperiment::colData(thresholded_data())), multiple = FALSE), 
                 selectInput("div_plot_over", label = "Plot Over", choices = colnames(SummarizedExperiment::colData(thresholded_data())), multiple = FALSE), 
                 textInput("div_your_title", "Title", value = ""),
                 br(),
                 numericInput("div_line_size", "Line size", value = 2),
                 numericInput("div_point_size", "Point size", value = 5),
                 numericInput("div_text_size", "Text size", value = 16)
               )
        ),
        
        column(8,
               plotOutput('view_div_plot', height = 500)
        )
      )
      
    })
    

    #======================================================================================================
    # Ridge Plot TAB
    output$RidgePlot <- renderUI({
      
      if (is.null(thresholded_data()))
        return()
      
      RidgePlotInput <- function(){
        print(barcodetrackR::ridge_plot(your_SE = thresholded_data(),
                                        split_bias_on = input$ridge_split_bias_on,
                                        bias_1 = input$ridge_bias1,
                                        bias_2 = input$ridge_bias2,
                                        split_bias_over = input$ridge_split_bias_over,
                                        remove_unique = input$ridge_rem_unique,
                                        weighted = input$ridge_weighted,
                                        text_size = input$ridge_text_size,
                                        add_dots = input$ridge_add_dots
        ))
      }
      
      
      output$view_ridge_plot <- renderPlot({
        RidgePlotInput()
      })
      
      # Helper function to get unique values of the split_bias_on parameter
      observeEvent(input$ridge_split_bias_on, {updateSelectizeInput(session,
                                                             inputId = 'ridge_bias1',
                                                             choices = unique(SummarizedExperiment::colData(thresholded_data())[,input$ridge_split_bias_on]))})
      
      observeEvent(input$ridge_split_bias_on, {updateSelectizeInput(session,
                                                             inputId = 'ridge_bias2',
                                                             choices = unique(SummarizedExperiment::colData(thresholded_data())[,input$ridge_split_bias_on]))})
      
      
      fluidRow(
        column(3,
               wellPanel(
                 selectInput("ridge_split_bias_on", label = "Compare Samples by", choices = colnames(SummarizedExperiment::colData(thresholded_data())), multiple = FALSE), 
                 selectInput("ridge_bias1", label = "Selection 1", choices = c(""), multiple = FALSE), 
                 selectInput("ridge_bias2", label = "Selection 2", choices = c(""), multiple = FALSE), 
                 selectInput("ridge_split_bias_over", label = "Plot Over", choices = colnames(SummarizedExperiment::colData(thresholded_data())), multiple = FALSE), 
                 strong("Options"),
                 checkboxInput("ridge_weighted", "Weighted", value = FALSE),
                 checkboxInput("ridge_rem_unique", "Remove unique barcodes", value = FALSE),
                 checkboxInput("ridge_add_dots", "Add dots", value = TRUE),
                 numericInput("ridge_text_size", "Text size", value = 16)
               )
        ),
        
        column(8,
               plotOutput('view_ridge_plot', height = 500)
        )
      )

      
    })
    
    
    #======================================================================================================
    # Chord Diagram TAB
    output$ChordDiagram <- renderUI({
      
      if (is.null(thresholded_data()))
        return()
      
      ChordDiagramInput <- function(){
        print(barcodetrackR::circos_plot(your_SE = circos_data(),
                                         weighted = input$circos_weighted,
                                         your_title = input$circos_title,
                                         plot_label = input$circos_plot_label,
                                         alpha = input$circos_alpha
        ))
      }
      
      circos_data <- reactive({
        se <- thresholded_data()
        se <- se[,se$SAMPLENAME %in% input$circos_Samples] # subset samples
        se$SAMPLENAME <- factor(se$SAMPLENAME, levels = input$circos_Samples)
        se <- se[,order(se$SAMPLENAME)]
        return(se)
      })
      
      output$view_circos_plot <- renderPlot({
        ChordDiagramInput()
      })
      
      fluidRow(
        column(3,
               wellPanel(
                 selectizeInput("circos_Samples", "Samples", choices = as.vector(unique(thresholded_data()$SAMPLENAME)), multiple = TRUE),
                 strong("Options"),
                 checkboxInput("circos_weighted", "Weighted", value = FALSE),
                 selectInput("circos_plot_label", "Plot labels",choices = colnames(SummarizedExperiment::colData(thresholded_data())), multiple = FALSE),
                 textInput("circos_title", "Title", value = ""),
                 numericInput("circos_alpha", "Alpha", value = 1)
               )
        ),
        
        column(8,
               plotOutput('view_circos_plot', height = 400)
        )
      )
      
    })
    
  
    
    # #======================================================================================================
    # #SCATTER PLOT TAB
    # 
    # output$ScatterPlot <- renderUI({
    #   if(is.null(thresholded_data()))
    #     return()
    # 
    #   scatterInput <- function(){
    #     normdata <- 100*prop.table(as.matrix(scatter_data()), margin = 2)
    #     plot(normdata, col = 'black', ylim = c(0,input$scatter_Ylim), xlim = c(0,input$scatter_Xlim), cex = input$scatter_Dotsize, pch = 20)
    #     legend("topright", legend = paste("Correlation (R) =", cor(normdata[,1], normdata[,2])))
    #   }
    # 
    #   output$viewscatter <- renderPlot({
    #     scatterInput()
    #     height = 700
    #   })
    # 
    #   scatter_data <- reactive({
    #     scf <- thresholded_data()
    # 
    #     scf <- scf[scf$GIVENNAME %in% c(input$scatter_Sample1, input$scatter_Sample2),] #subset samples
    #     newcolnames <- scf$GIVENNAME
    # 
    #     scf$GIVENNAME <- NULL
    #     scf$EXPERIMENT <- NULL
    #     scf$CELLTYPE <- NULL
    #     scf$MONTH <- NULL
    #     scf$LOCATION <- NULL
    #     scf$MISC <- NULL
    # 
    #     scf <- data.frame(t(scf))
    #     colnames(scf) <- newcolnames
    #     return(scf)
    # 
    #   })
    # 
    #   fluidRow(
    #     column(3,
    #            wellPanel(
    #              selectInput("scatter_Sample1", label = "1. First Sample to Use", choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = FALSE),
    #              selectInput("scatter_Sample2", label = "2. Second Sample to Use", choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = FALSE),
    #              numericInput("scatter_Threshold", "3. Select Threshold", value = 0, step = 100, min = 0),
    #              numericInput("scatter_Ylim", "4. Y Limit: ", value = 100, step = 0.1, min = 0),
    #              numericInput("scatter_Xlim", "5. X Limit: ", value = 100, step = 0.1, min = 0),
    #              numericInput("scatter_Dotsize","6. Dot Size: ", value = 5)
    # 
    #            )),
    # 
    # 
    # 
    #     column(6,
    #            plotOutput('viewscatter', height = 700)
    #     )
    # 
    #   )
    # 
    # 
    # })
    # 
    
    # 
    # 
    # # #======================================================================================================
    # #Binary Heatmap TAB
    # 
    # output$BinaryHeatmap <- renderUI({
    # 
    # 
    #   if (is.null(thresholded_data()))
    #     return()
    # 
    # 
    #   BinaryheatmapInput <- function(){
    #     print(barcodetrackR::barcode_binary_heatmap(your_data = Binaryheatmap_data(),
    #                                                 label_size = input$Binaryheatmap_labels,
    #                                                 percent_threshold = input$Binaryheatmap_threshold
    #     ))
    #   }
    # 
    #   output$viewBinaryheatmap <- renderPlot({
    #     BinaryheatmapInput()
    #     height = 900
    #   })
    # 
    # 
    #   Binaryheatmap_data <- reactive({
    #     BBdf <- thresholded_data()
    #     BBdf <- BBdf[BBdf$GIVENNAME %in% input$Binaryheatmap_samples,] #subset samples
    #     BBdf$GIVENNAME <- factor(BBdf$GIVENNAME, levels = input$Binaryheatmap_samples)
    #     BBdf <- BBdf[order(BBdf$GIVENNAME),]
    #     newcolnames <- BBdf$GIVENNAME
    #     BBdf$GIVENNAME <- NULL
    #     BBdf$EXPERIMENT <- NULL
    #     BBdf$CELLTYPE <- NULL
    #     BBdf$MONTH <- NULL
    #     BBdf$LOCATION <- NULL
    #     BBdf$MISC <- NULL
    #     BBdf <- data.frame(t(BBdf))
    #     colnames(BBdf) <- newcolnames
    #     return(BBdf)
    #   })
    # 
    #   fluidRow(
    #     column(3,
    #            wellPanel(
    #              selectizeInput("Binaryheatmap_samples", label = "1. Which Samples to Use (in order)",
    #                             choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = TRUE),
    #              numericInput("Binaryheatmap_threshold", "2. Set Threshold", value = 0),
    #              numericInput("Binaryheatmap_labels", "3. Set Column Label Size", value = 15)
    #            )
    #     ),
    #     column(8,
    #            plotOutput('viewBinaryheatmap', height = 700)
    #     )
    # 
    #   )
    # 
    # })
    # 
    # #======================================================================================================
    # #TERNPLOT TAB
    # 
    # output$TernPlot <- renderUI({
    #   if(is.null(thresholded_data()))
    #     return()
    # 
    #   ternplotInput <- function(){
    #     print(barcodetrackR::ternary_plot(ternplot_data(),
    #                                       show_arrows = input$ternplot_Showarrows,
    #                                       show_ticks = input$ternplot_Showticks,
    #                                       density_mode = input$ternplot_Density))
    # 
    #   }
    # 
    #   output$viewternplot<- renderPlot({
    #     ternplotInput()
    #     height = 700
    #   })
    # 
    #   ternplot_data <- reactive({
    #     terndf <- thresholded_data()
    # 
    #     terndf <- terndf[terndf$GIVENNAME %in% input$ternplot_Samples,] #subset samples
    #     terndf$GIVENNAME <- factor(terndf$GIVENNAME, levels = input$ternplot_Samples)
    #     terndf <- terndf[order(terndf$GIVENNAME),]
    #     newcolnames <- terndf$GIVENNAME
    # 
    #     terndf$GIVENNAME <- NULL
    #     terndf$EXPERIMENT <- NULL
    #     terndf$CELLTYPE <- NULL
    #     terndf$MONTH <- NULL
    #     terndf$LOCATION <- NULL
    #     terndf$MISC <- NULL
    # 
    #     terndf <- data.frame(t(terndf))
    #     colnames(terndf) <- newcolnames
    #     return(terndf)
    # 
    #   })
    # 
    #   fluidRow(
    #     column(3,
    #            wellPanel(
    #              selectInput("ternplot_Samples", label = "1. Which Samples to Use (3)",
    #                          choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = TRUE),
    #              strong("2. Options"),
    #              checkboxInput("ternplot_Density", label = "Density Mode", value = FALSE),
    #              checkboxInput("ternplot_Showarrows", label = "Show Arrows", value = TRUE),
    #              checkboxInput("ternplot_Showticks", label = "Show Tick Marks", value = FALSE)
    #            )),
    # 
    # 
    # 
    #     column(9,
    #            plotOutput('viewternplot', height = 700)
    #     )
    # 
    #   )
    # 
    # 
    # 
    # 
    # 
    # 
    # })
    # 
    # 
    # #======================================================================================================
    # #RANKABUNDANCE TAB
    # 
    # output$RankAbundance <- renderUI({
    #   if(is.null(thresholded_data()))
    #     return()
    # 
    #   rankabundanceInput <- function(){
    #     print(barcodetrackR::rank_abundance_plot(rankabundance_data(),
    #                                              dot_size = input$rankabundance_Dotsize,
    #                                              text_size = input$rankabundance_Textsize))
    # 
    #   }
    # 
    #   output$viewrankabundance<- renderPlot({
    #     rankabundanceInput()
    #     height = 700
    #   })
    # 
    #   rankabundance_data <- reactive({
    #     rankabdf <- thresholded_data()
    #     rankabdf <- rankabdf[rankabdf$GIVENNAME %in% input$rankabundance_Samples,] #subset samples
    #     rankabdf$GIVENNAME <- factor(rankabdf$GIVENNAME, levels = input$rankabundance_Samples)
    #     rankabdf <- rankabdf[order(rankabdf$GIVENNAME),]
    #     newcolnames <- rankabdf$GIVENNAME
    # 
    #     rankabdf$GIVENNAME <- NULL
    #     rankabdf$EXPERIMENT <- NULL
    #     rankabdf$CELLTYPE <- NULL
    #     rankabdf$MONTH <- NULL
    #     rankabdf$LOCATION <- NULL
    #     rankabdf$MISC <- NULL
    # 
    #     rankabdf <- data.frame(t(rankabdf))
    #     colnames(rankabdf) <- newcolnames
    #     return(rankabdf)
    # 
    #   })
    # 
    #   fluidRow(
    #     column(3,
    #            wellPanel(
    #              selectInput("rankabundance_Samples", label = "1. Which Samples to Use",
    #                          choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = TRUE),
    #              numericInput("rankabundance_Dotsize", "2. Enter Dot Size: ", value = 3),
    #              numericInput("rankabundance_Textsize", "3. Enter Text Size: ", value = 15)
    #            )),
    # 
    # 
    # 
    #     column(9,
    #            plotOutput('viewrankabundance', height = 700)
    #     )
    # 
    #   )
    # 
    # 
    # 
    # 
    # 
    # 
    # })
    # 
    # #======================================================================================================
    # #CLONALBIAS TAB
    # 
    # output$ClonalBias <- renderUI({
    #   if(is.null(thresholded_data()))
    #     return()
    # 
    #   clonalbiasInput <- function(){
    #     if(input$clonalbias_Type == "Dots"){
    #       print(barcodetrackR::dot_bias(clonalbias_data(), text_size = input$clonalbias_Textsize))
    #     } else if(input$clonalbias_Type == "Bars"){
    #       print(barcodetrackR::clonal_bias(clonalbias_data(), text_size = input$clonalbias_Textsize))
    #     }
    #   }
    # 
    #   output$viewclonalbias<- renderPlot({
    #     clonalbiasInput()
    #     height = 700
    #   })
    # 
    #   clonalbias_data <- reactive({
    #     cb_df <- thresholded_data()
    #     cb_df <- cb_df[cb_df$GIVENNAME %in% input$clonalbias_Samples,] #subset samples
    #     cb_df$GIVENNAME <- factor(cb_df$GIVENNAME, levels = input$clonalbias_Samples)
    #     cb_df <- cb_df[order(cb_df$GIVENNAME),]
    #     newcolnames <- cb_df$GIVENNAME
    #     cb_df$GIVENNAME <- NULL
    #     cb_df$EXPERIMENT <- NULL
    #     cb_df$CELLTYPE <- NULL
    #     cb_df$MONTH <- NULL
    #     cb_df$LOCATION <- NULL
    #     cb_df$MISC <- NULL
    #     cb_df <- data.frame(t(cb_df))
    #     colnames(cb_df) <- newcolnames
    #     return(cb_df)
    #   })
    #   fluidRow(
    #     column(3,
    #            wellPanel(
    #              selectInput("clonalbias_Samples", label = "1. Which Samples to Use",
    #                          choices = as.vector(unique(thresholded_data()$GIVENNAME)), multiple = TRUE),
    #              numericInput("clonalbias_Textsize", "2. Enter Text Size: ", value = 15),
    #              selectInput("clonalbias_Type", label = "3. Pick type of plot",
    #                          choices = c("Dots", "Bars"))
    #            )),
    #     column(9,
    #            plotOutput('viewclonalbias', height = 700)
    #     )
    #   )
    # })






  }
)





