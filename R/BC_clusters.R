#'@title BC_clusters (Barcode Clusters)
#'
#'@description Using the same methodlogy as BCheatmap, takes the top 'n' rows from each column, clusters them, and plots
#'their average values over each sample (over time).
#'
#'@param your_data A data frame. Usually individual barcodes in rows and samples in columns.
#'@param names Names for each column of the data frame.
#'@param n_clones The top 'n' clones to plot.
#'@param your_title The title for the plot.
#'@param log_transform Logical. Log transform data before clustering and plotting.
#'@param log_choice Data is log transformed with this log.
#'@param variable_log_min Logical. TRUE to make log(0) = minimum nonzero - 1. FALSE to make log(0) = log(100/4000000) - 1.
#'@param distance_method Character. Use summary(proxy::pr_DB) to see all options.
#'@param minkowski_power The power of the Minkowski distance (if minkowski is used).
#'@param hclust_linkage Character. One of one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#'@param clusters How many clusters to cut hierarchical tree into for display when row_order is "hierarchical".
#'@param colors Which colors to plot (default to all, can choose using BCheatmap, unless using "BCheatmap")
#'@return Displays a line plot in the current plot window.
#'@examples
#'BCheatmap(your_data = zh33, names = colnames(zh33), n_clones = 10,
#'       your_title = "First Time Point", grid = TRUE, columnLabels = 3)
#'BCheatmap(your_data = zh33, n_clones = 10, printtable = TRUE)
#'@export

BC_clusters <- function(your_data,
                        names = colnames(your_data),
                        n_clones = 10,
                        months = 1:ncol(your_data),
                        your_title = "",
                        log_transform = TRUE,
                        log_choice = exp(1),
                        variable_log_min = TRUE,
                        distance_method = "Euclidean",
                        minkowski_power = 1,
                        hclust_linkage = "complete",
                        clusters = 3,
                        colors = "",
                        percent_scale = FALSE) {
  
  if(clusters <= 0 || colors == ""){
    stop("Must specify at least 1 cluster and 1 color to be plotted")
  }
  
  colnames(your_data) <- months
  
  #scales all data to be a percentage of reads instead of number of reads and keeps copy of raw read number
  your_data <- as.data.frame(prop.table(as.matrix(your_data),2))
  your_data[is.na(your_data)] <- 0
  if (any(colSums(your_data) != 1)){
    stop("One of your columns contained no data")
  }
  
  #creates data frame that shows rank of original your_data
  your_data_ranked <- apply(-your_data, 2, rank, ties.method = "min", na.last = "keep")
  
  #subsets those barcodes that have at least one top N clone
  top_clones_choices <- apply(your_data_ranked, 1, function(x){any(x<=n_clones, na.rm = TRUE)})
  your_data <- your_data[top_clones_choices,]
  your_data[is.na(your_data)] <- 0
  
  #takes log of data
  your_data_logged <- custom_log(your_data, log_choice, variable_log_min)
  
  if(distance_method == "Minkowski"){
    hclustering <-hclust(proxy::dist(your_data_logged, method = distance_method, p = minkowski_power), method = hclust_linkage)
  } else {
    hclustering <-hclust(proxy::dist(your_data_logged, method = distance_method), method = hclust_linkage)
  }
  e <- rev(hclustering$order)
  
  myPalette <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
                 "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
                 "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D",
                 "#8A7C64", "#599861")
  names(myPalette) <- c("Baby Blue", "Orange", "Neon Green", "Fuschia", "Ew Green", "Faint Red",
                        "Yucky Yellow", "Royal Blue", "Dark Purple", "Radioactive Yellow", "Purple Black",
                        "Faded Green", "Tan", "Oscar the Grouch Green", "Brown", "Old Bubble Gum", "Vivid Pink",
                        "80s Green", "Dark Brown", "Beach Blue", "Crimson Red", "Triceratops Purple", "Cool Blue",
                        "Gold", "Grandpa Beige", "Pine Green")
  ourZero = (log(100/4000000) - 1)
  
  cuts <- cutree(hclustering, k = clusters)
  cluster_colors = myPalette[1:clusters][cuts]
  std_error <- function(x) sd(x)/sqrt(length(x))
  final_data <- if(log_transform) your_data_logged else your_data
  melted_data <- reshape2::melt(cbind(cluster_colors[e], (final_data[e,])), id.vars = c(1))
  colnames(melted_data) <- c("COLOR", "SAMPLE", "PROP")
  if(percent_scale){
    melted_data$PROP <- ifelse(melted_data$PROP == ourZero, 0, exp(melted_data$PROP))
  }
  meltySD <- aggregate(melted_data$PROP, by = list(melted_data$COLOR, melted_data$SAMPLE), FUN = std_error)
  melted_data <- aggregate(melted_data$PROP, by = list(melted_data$COLOR, melted_data$SAMPLE), FUN = mean)
  #print(aggregate(melted_data$PROP, by = list(melted_data$COLOR, melted_data$SAMPLE), FUN = qt))
  melted_data$SD <- meltySD$x
  colnames(melted_data) <- c("COLOR", "SAMPLE", "PROP", "SD")
  colors <- myPalette[colors]
  print(myPalette[(unique(melted_data$COLOR))])
  melted_data <- melted_data[melted_data$COLOR %in% colors,]
  myBreaks = c(ourZero, pretty(min(melted_data$PROP):(ceiling(max(melted_data$PROP)))))
  if(percent_scale){
    myBreaks = pretty(1:ceiling(max(melted_data$PROP)))
    print(myBreaks)
  }
  #myBreaks <- myBreaks[myBreaks >= ourZero]
  formatBack <- function(x){paste0(x*100, " %")}
  melted_data$SAMPLE <- as.numeric(as.character(melted_data$SAMPLE))
  ggplot2::ggplot(melted_data, ggplot2::aes(x = SAMPLE, y = PROP, group = COLOR))+
    ggplot2::geom_point(ggplot2::aes(color = COLOR), size = 4)+
    #ggplot2::scale_y_continuous(breaks = myBreaks, labels = formatBack)+
    #ggplot2::scale_y_continuous(breaks = myBreaks, labels = formatBack)+
    ggplot2::scale_y_continuous(breaks = c(0,0.02, 0.04, 0.06, 0.08, 0.1), labels = formatBack)+
    ggplot2::scale_x_continuous(breaks = melted_data$SAMPLE)+
    ggplot2::scale_color_identity()+
    ggplot2::geom_line(ggplot2::aes(color = COLOR), size = 3)+
    ggplot2::geom_ribbon(ggplot2::aes(ymin = PROP - SD, ymax = PROP + SD, fill = COLOR, color = COLOR), alpha = 0.3, show.legend = FALSE)+
    ggplot2::scale_fill_identity()+
    ggplot2::ggtitle(paste0('\n',your_title, '\n'))+
    ggplot2::xlab("\nMonth\n")+
    ggplot2::ylab("\n%\n")+
    ggplot2::coord_cartesian(ylim = c(0,.04))+
    ggplot2::theme(text = ggplot2::element_text(size=30),
                   panel.background = ggplot2::element_rect(fill = "white", colour = "black"),
                   panel.grid.major = ggplot2::element_line(colour = "lightgrey"),
                   axis.text.x = ggplot2::element_text(angle = 90))
}

custom_log <- function(x, log_choice, vary_log_min){
  x <- log(x, log_choice)
  #sets the -inf values to be the minimum of the possible value in the data frame -1
  if(vary_log_min) {
    x[x == -Inf] <- (min(x[x > -Inf]) - 1)
  } else {
    x[x == -Inf] <- (log(100/4000000) - 1)
  }
  return(x)
}

