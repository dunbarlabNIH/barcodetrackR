#' Barcode Circos plot
#'
#' Creates a chord diagram showing each cell type as a region around a circle and shared clones between these cell types as links between the regions.
#'
#'@param your_SE A Summarized Experiment object.
#'@param weighted weighted = F which is default will make links based on # of shared clones. Weighted = T will make links based on their proportion.
#'@param plot_label Name of colData variable to use as labels for regions. Defaults to SAMPLENAME
#'@param alpha Transparency of links. Default = 1 is opaque. 0 is completely transluscent
#'@param your_title The title for the plot.
#'@param text_size Size of region labels
#'
#'@return Displays a chord diagram in the current plot window.
#'
#'@import viridis
#'@importFrom dplyr group_by_all
#'
#'@export
#'
#'@examples
#'circos_plot(your_SE = ZH33_SE,  plot_label = 'Cell_type')
#'
circos_plot <- function(your_SE,
                        weighted = FALSE,
                        plot_label = "SAMPLENAME",
                        alpha = 1,
                        your_title = NULL,
                        text_size = 12) {

# Load data, remove data that is zero in all timepoints
your_data <- SummarizedExperiment::assays(your_SE)$counts
meta_data <- SummarizedExperiment::colData(your_SE)
your_data <- as.matrix(your_data[rowSums(your_data) > 0,])

#get labels for heatmap
# plot_label <- plot_label %||% 'SAMPLENAME'
colnames(your_data) <- meta_data[,plot_label]

# Create binary matrix of data
temp_binary <- as.matrix((your_data > 0)+0)
# Create temp of proportions
temp_prop <- prop.table(your_data,2)*100

# Sort temp prop by binary temp
temp_prop <- as.matrix(temp_prop[do.call(order,-as.data.frame(temp_binary)),])
# Sort binary temp
temp_binary <- as.matrix(temp_binary[do.call(order,-as.data.frame(temp_binary)),])

# Remove solo clones, they are just the space left at the end
temp_prop <- temp_prop[rowSums(temp_binary)>1,]
temp_binary <- temp_binary[rowSums(temp_binary)>1,]

# Count the number of occurences of the unique combinations
unique_count <- as.data.frame(temp_binary) %>% dplyr::group_by_all %>% count
# Sort decreasing
unique_count <- unique_count[do.call(order,-unique_count),]

# Get the proportions of each barcode matching each unique combination
unique_prop <- unique_count
unique_prop$n <- NULL
count_vec <- unique_count$n
my_counter <- 1
for (i in 1:nrow(unique_count)){
  my_start <- my_counter
  my_end <- my_counter + count_vec[i] - 1
  if (my_start == my_end){
    unique_prop[i,] <- as.list(temp_prop[my_start:my_end,])
  } else {
    unique_prop[i,] <- as.list(colSums(temp_prop[my_start:my_end,]))
  }
  my_counter <- my_counter + as.numeric(unique_count$n[i])
}

# Generate counting index for each cell type
count_index <- rep(0,length(colnames(temp_binary)))
names(count_index) <- colnames(temp_binary)

prop_count_index <- rep(0,length(colnames(temp_prop)))
names(prop_count_index) <- colnames(temp_binary)

# Set up color pallete
my_cols = viridis::viridis(nrow(unique_count))


if (weighted == FALSE){
  # Initialize circos plot
  par(mar = c(0.5, 0.5, 1, 0.5))
  circlize::circos.par(points.overflow.warning = FALSE)
  # circos.par(cell.padding = c(0, 0, 0, 0))

  # Make x limits matrix
  xlims <- matrix(data = 0, nrow = length(colnames(your_data)), ncol = 2) # Just a placeholder
  for (m in 1:length(colnames(your_data))){
    xlims[m,2] <- colSums(your_data != 0)[m]
  }
  circlize::circos.initialize(factors = factor(colnames(your_data), levels = colnames(your_data)), xlim = xlims)

  # Create outer tracks of circos plot
  circlize::circos.track(factors = factor(colnames(your_data), levels = colnames(your_data)), ylim = c(0, 1), bg.col = "grey",
               bg.border = NA, track.height = 0.1)
  # Add labels
  circlize::circos.trackText(x = xlims[,2]/2, y = rep(0.5,length(colnames(your_data))),
                             factors = factor(colnames(your_data),levels = colnames(your_data)),
                             labels = factor(colnames(your_data), levels = colnames(your_data)),
                             niceFacing = T, cex = text_size/12)

  # Loop through rows of unique_count
  for (i in 1:nrow(unique_count)){
    num_cells <- sum(unique_count[i,1:length(colnames(temp_binary))])
    num_links <- num_cells*(num_cells-1)/2

    cell_list <- colnames(temp_binary)[which(unique_count[i,1:length(colnames(temp_binary))]>0)]
    comb_mat <- combn(cell_list,2)
    # Loop through number of links that must be drawn
    for (j in 1:num_links){
      # Draw links
      cell.1 <- comb_mat[1,j]
      cell.2 <- comb_mat[2,j]
      circlize::circos.link(cell.1, c(count_index[cell.1],count_index[cell.1] + as.numeric(unique_count[i,"n"])),
                            cell.2, c(count_index[cell.2],count_index[cell.2] + as.numeric(unique_count[i,"n"])),
                            col = adjustcolor(my_cols[i],alpha.f = alpha))
    }
    # Update indices
    for (k in 1:num_cells){
      count_index[cell_list[k]] <- count_index[cell_list[k]] + as.numeric(unique_count[i,"n"])
    }
  }
  title(your_title, adj = 0)
  # circos.clear()
}

else if (weighted == TRUE){
  # Initialize circos plot
  par(mar = c(0.5, 0.5, 1, 0.5))
  circlize::circos.par(points.overflow.warning = FALSE)

  # Make x limits matrix
  xlims <- t(replicate(length(colnames(your_data)), c(0,100)))
  circlize::circos.initialize(factors = factor(colnames(your_data), levels = colnames(your_data)), xlim = xlims)

  # Create outer tracks of circos plot
  circlize::circos.track(factors = factor(colnames(your_data), levels = colnames(your_data)), ylim = c(0, 1), bg.col = "grey",
                         bg.border = NA, track.height = 0.1)
  # Add labels
  circlize::circos.trackText(x = xlims[,2]/2, y = rep(0.5,length(colnames(your_data))),
                             factors = factor(colnames(your_data),levels = colnames(your_data)),
                             labels = factor(colnames(your_data), levels = colnames(your_data)),
                             niceFacing = T, cex = text_size/12)

  # Loop through rows of unique_count
  for (i in 1:nrow(unique_prop)){
    num_cells <- sum(unique_count[i,1:length(colnames(temp_prop))])
    num_links <- num_cells*(num_cells-1)/2

    cell_list <- colnames(temp_prop)[which(unique_prop[i,1:length(colnames(temp_prop))]>0)]
    comb_mat <- combn(cell_list,2)
    # Loop through number of links that must be drawn
    for (j in 1:num_links){
      # Draw links
      cell.1 <- comb_mat[1,j]
      cell.2 <- comb_mat[2,j]
      circlize::circos.link(cell.1, c(prop_count_index[cell.1],prop_count_index[cell.1] + as.numeric(unique_prop[i,cell.1])),
                            cell.2, c(prop_count_index[cell.2],prop_count_index[cell.2] + as.numeric(unique_prop[i,cell.2])),
                            col = adjustcolor(my_cols[i],alpha.f = alpha))
    }
    # Update indices
    for (k in 1:num_cells){
      prop_count_index[cell_list[k]] <- prop_count_index[cell_list[k]] + as.numeric(unique_prop[i,cell_list[k]])
    }
  }
  title(your_title, adj = 0)
  # circlize::circos.clear()
}
my_p  <-par()
circlize::circos.clear()

invisible(my_p)
}
