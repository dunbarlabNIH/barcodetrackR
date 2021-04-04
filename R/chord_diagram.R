#' Barcode Chord Diagram
#'
#' Creates a chord diagram showing each cell type (or other factor) as a region around a circle and shared clones between these cell types as links between the regions. The space around the regions which is not connected to a chord indicates clones unique to that sample, not shared with other samples.
#'
#' @param your_SE Summarized Experiment object containing clonal tracking data as created by the barcodetrackR `create_SE` function.
#' @param weighted Logical. weighted = F which is default will make links based on the number of shared clones between the factors. Weighted = TRUE will make the link width based on the clone's proportion in the samples.
#' @param plot_label Character. Name of colData variable to use as labels for regions. Defaults to SAMPLENAME
#' @param alpha Numeric. Transparency of links. Default = 1 is opaque. 0 is completely transluscent
#' @param your_title Character. The title for the plot.
#' @param text_size Numeric. Size of region labels
#' @param return_table Logical. If set to TRUE, in addition to plotting a chord diagram in the plot window, the function will return a dataframe of the shared clonality used to make the chord diagram. If Weighted is FALSE, the dataframe will contain a row for each set of clones and values of 1 or 0 indicating the samples which share that set of clones, and a freq column which is the number of clones in that set. If weighted is set to TRUE, each row will contain a set of clones and the data will show the proportion that set of clones comprises in each sample. The proportions of 0 indicate which samples do not share that set of clones.
#'
#' @return Displays a chord diagram in the current plot window depicting shared clonality between samples (regions) as chords or links between the regions. Or,
#'
#' @import viridis
#' @import graphics
#' @importFrom grDevices adjustcolor
#'
# #'@import dplyr
#' @rawNamespace import(dplyr, except = count)
#' @importFrom plyr count
#'
#' @export
#'
#' @examples
#' data(wu_subset)
#' chord_diagram(your_SE = wu_subset[, c(4, 8, 12)], plot_label = "celltype")
chord_diagram <- function(your_SE,
    weighted = FALSE,
    plot_label = "SAMPLENAME",
    alpha = 1,
    your_title = NULL,
    text_size = 12,
    return_table = FALSE) {

    # Load data, remove data that is zero in all timepoints
    your_data <- SummarizedExperiment::assays(your_SE)$counts
    meta_data <- SummarizedExperiment::colData(your_SE)
    your_data <- as.matrix(your_data[rowSums(your_data) > 0, ])

    # Error check
    if (!plot_label %in% colnames(meta_data)) {
        stop("Provided plot_label is not a column of the metadata.")
    }
    # get labels for heatmap
    colnames(your_data) <- meta_data[, plot_label]

    # Create binary matrix of data
    temp_binary <- as.matrix((your_data > 0) + 0)
    # Create temp of proportions
    temp_prop <- prop.table(your_data, 2) * 100

    # Sort temp prop by binary temp
    temp_prop <- as.matrix(temp_prop[do.call(order, -as.data.frame(temp_binary)), ])
    # Sort binary temp
    temp_binary <- as.matrix(temp_binary[do.call(order, -as.data.frame(temp_binary)), ])

    # Remove solo clones, they are just the space left at the end
    temp_prop <- temp_prop[rowSums(temp_binary) > 1, ]
    temp_binary <- temp_binary[rowSums(temp_binary) > 1, ]

    # Count the number of occurences of the unique combinations
    unique_count <- plyr::count(as.data.frame(temp_binary))
    # Sort decreasing
    unique_count <- unique_count[do.call(order, -unique_count), ]

    # Get the proportions of each barcode matching each unique combination
    unique_prop <- unique_count
    unique_prop$freq <- NULL
    count_vec <- unique_count$freq
    my_counter <- 1
    for (i in seq_len(nrow(unique_count))) {
        my_start <- my_counter
        my_end <- my_counter + count_vec[i] - 1
        if (my_start == my_end) {
            unique_prop[i, ] <- as.list(temp_prop[my_start:my_end, ])
        } else {
            unique_prop[i, ] <- as.list(colSums(temp_prop[my_start:my_end, ]))
        }
        my_counter <- my_counter + as.numeric(unique_count$freq[i])
    }

    # Generate counting index for each cell type
    count_index <- rep(0, length(colnames(temp_binary)))
    names(count_index) <- colnames(temp_binary)

    prop_count_index <- rep(0, length(colnames(temp_prop)))
    names(prop_count_index) <- colnames(temp_binary)

    # Set up color pallete
    my_cols <- viridis::viridis(nrow(unique_count))


    if (!weighted) {
        if (return_table) {
            return(unique_count)
        }

        # Initialize circos plot
        par(mar = c(0.5, 0.5, 1, 0.5))
        circlize::circos.par(points.overflow.warning = FALSE)

        # Make x limits matrix
        xlims <- matrix(data = 0, nrow = length(colnames(your_data)), ncol = 2) # Just a placeholder
        for (m in seq_along(colnames(your_data))) {
            xlims[m, 2] <- colSums(your_data != 0)[m]
        }
        circlize::circos.initialize(factors = factor(colnames(your_data), levels = colnames(your_data)), xlim = xlims)

        # Create outer tracks of circos plot
        circlize::circos.track(
            factors = factor(colnames(your_data), levels = colnames(your_data)), ylim = c(0, 1), bg.col = "grey",
            bg.border = NA, track.height = 0.1
        )
        # Add labels
        circlize::circos.trackText(
            x = xlims[, 2] / 2, y = rep(0.5, length(colnames(your_data))),
            factors = factor(colnames(your_data), levels = colnames(your_data)),
            labels = factor(colnames(your_data), levels = colnames(your_data)),
            niceFacing = TRUE, cex = text_size / 12
        )

        # Loop through rows of unique_count
        for (i in seq_len(nrow(unique_count))) {
            num_cells <- sum(unique_count[i, seq_along(colnames(temp_binary))])
            num_links <- num_cells * (num_cells - 1) / 2

            cell_list <- colnames(temp_binary)[which(unique_count[i, seq_along(colnames(temp_binary))] > 0)]
            comb_mat <- combn(cell_list, 2)
            # Loop through number of links that must be drawn
            for (j in seq_len(num_links)) {
                # Draw links
                cell.1 <- comb_mat[1, j]
                cell.2 <- comb_mat[2, j]
                circlize::circos.link(cell.1, c(count_index[cell.1], count_index[cell.1] + as.numeric(unique_count[i, "freq"])),
                    cell.2, c(count_index[cell.2], count_index[cell.2] + as.numeric(unique_count[i, "freq"])),
                    col = adjustcolor(my_cols[i], alpha.f = alpha)
                )
            }
            # Update indices
            for (k in seq_len(num_cells)) {
                count_index[cell_list[k]] <- count_index[cell_list[k]] + as.numeric(unique_count[i, "freq"])
            }
        }
        title(your_title, adj = 0)
    }

    else if (weighted) {
        if (return_table) {
            return(unique_prop)
        }

        # Initialize circos plot
        par(mar = c(0.5, 0.5, 1, 0.5))
        circlize::circos.par(points.overflow.warning = FALSE)

        # Make x limits matrix
        xlims <- t(replicate(length(colnames(your_data)), c(0, 100)))
        circlize::circos.initialize(factors = factor(colnames(your_data), levels = colnames(your_data)), xlim = xlims)

        # Create outer tracks of circos plot
        circlize::circos.track(
            factors = factor(colnames(your_data), levels = colnames(your_data)), ylim = c(0, 1), bg.col = "grey",
            bg.border = NA, track.height = 0.1
        )
        # Add labels
        circlize::circos.trackText(
            x = xlims[, 2] / 2, y = rep(0.5, length(colnames(your_data))),
            factors = factor(colnames(your_data), levels = colnames(your_data)),
            labels = factor(colnames(your_data), levels = colnames(your_data)),
            niceFacing = TRUE, cex = text_size / 12
        )

        # Loop through rows of unique_count
        for (i in seq_len(nrow(unique_prop))) {
            num_cells <- sum(unique_count[i, seq_along(colnames(temp_prop))])
            num_links <- num_cells * (num_cells - 1) / 2

            cell_list <- colnames(temp_prop)[which(unique_prop[i, seq_along(colnames(temp_prop))] > 0)]
            comb_mat <- combn(cell_list, 2)
            # Loop through number of links that must be drawn
            for (j in seq_len(num_links)) {
                # Draw links
                cell.1 <- comb_mat[1, j]
                cell.2 <- comb_mat[2, j]
                circlize::circos.link(cell.1, c(prop_count_index[cell.1], prop_count_index[cell.1] + as.numeric(unique_prop[i, cell.1])),
                    cell.2, c(prop_count_index[cell.2], prop_count_index[cell.2] + as.numeric(unique_prop[i, cell.2])),
                    col = adjustcolor(my_cols[i], alpha.f = alpha)
                )
            }
            # Update indices
            for (k in seq_len(num_cells)) {
                prop_count_index[cell_list[k]] <- prop_count_index[cell_list[k]] + as.numeric(unique_prop[i, cell_list[k]])
            }
        }
        title(your_title, adj = 0)
    }
    my_p <- par()
    circlize::circos.clear()

    invisible(my_p)
}
