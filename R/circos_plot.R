#' Barcode Circos plot
#'
#' Creates a circos plot showing each cell type as a region around a circle and shared clones between these cell types as links between the regions.
#'
#'@param your_SE A Summarized Experiment object.
#'@param plot_label Name of colData variable to use as labels for regions. Defaults to SAMPLENAME
#'
#'@return Displays a circos plot in the current plot window.
#'
#'@export
#'
#'@examples
#'barcode_ggheatmap_2(your_SE = ZH33_SE,  plot_label = 'Cell_type')
#'
circos_plot <- function(your_SE,
                        plot_label = "SAMPLENAME",
                        your_title = "",
                        color_scale = NULL) {
  
# Remove data that is zero in all timepoints
your_data <- SummarizedExperiment::assays(your_SE)$counts
meta_data <- SummarizedExperiment::colData(your_SE)
your_data <- your_data[rowSums(your_data) > 0,]

#get labels for heatmap
# plot_label <- plot_label %||% 'SAMPLENAME'
colnames(your_data) <- meta_data[,plot_label]

# Initialize circos plot
circos.par(points.overflow.warning = FALSE)

# Make x limits matrix
xlims <- matrix(data = 0, nrow = length(colnames(your_data)), ncol = 2) # Just a placeholder
for (m in 1:length(colnames(your_data))){
  xlims[m,2] <- colSums(your_data != 0)[m]
}
circos.initialize(factors = factor(colnames(your_data), levels = colnames(your_data)), xlim = xlims)

# Create outer tracks of circos plot
circos.track(factors = factor(colnames(your_data), levels = colnames(your_data)), ylim = c(0, 1), bg.col = "grey",
             bg.border = NA, track.height = 0.1)
# Add labels
circos.trackText(x = xlims[,2]/2, y = rep(0.5,length(colnames(your_data))),
                 factors = factor(colnames(your_data),levels = colnames(your_data)),
                 labels = factor(colnames(your_data), levels = colnames(your_data)),
                 niceFacing = T)

# Create binary matrix of data
temp_binary <- as.matrix((your_data > 0)+0)
# Sort binary temp
temp_binary <- as.matrix(temp_binary[do.call(order,-as.data.frame(temp_binary)),])
# Remove solo clones, they are just the space left at the end
temp_binary <- temp_binary[rowSums(temp_binary)>1,]
# Create matrix of the unique combinations
binary_unique <- unique(temp_binary)
# Count the number of occurences of the unique combinations
unique_count <- as.data.frame(temp_binary) %>% group_by_all %>% count
# Sort decreasing
unique_count <- unique_count[do.call(order,-unique_count),]

# Generate counting index for each cell type
count_index <- rep(0,length(colnames(temp_binary)))
names(count_index) <- colnames(temp_binary)

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
    circos.link(cell.1, c(count_index[cell.1],count_index[cell.1] + unique_count[i,"freq"]), 
                cell.2, c(count_index[cell.2],count_index[cell.2] + unique_count[i,"freq"]), 
                col = i)
  }
  # Update indices
  for (k in 1:num_cells){
    count_index[cell_list[k]] <- count_index[cell_list[k]] + unique_count[i,"freq"]
  }
}

circos.clear()
}
