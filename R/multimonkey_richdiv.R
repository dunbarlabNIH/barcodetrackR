multimonkey_richness_plot <- function(outfile_list,
                                      outfile_names,
                                      outfile_samples_list,
                                      outfile_months_list,
                                      celltypes_list,
                                      threshold_list,
                                      combine = FALSE,
                                      richness_type = "unique",
                                      diversity_index = "shannon",
                                      point_size = 5, line_size = 3,
                                      y_lower = 0, y_upper = 9.5){

  if(combine){
    for(i in 1:length(outfile_list)){

      outfile_list[[i]] <- outfile_list[[i]][,outfile_samples_list[[i]], drop = FALSE]
      outfile_list[[i]][outfile_list[[i]] < threshold_list[[i]]] <- 0

      splitlist <- split(as.data.frame(t(outfile_list[[i]])), rep(1:length(outfile_months_list[[i]]), each = length(celltypes_list[[i]])))
      outfile_list[[i]] <- sapply(splitlist, colSums)
      colnames(outfile_list[[i]]) <- outfile_months_list[[i]]

      if(richness_type == "diversity"){
        stop("Cannot use combine == TRUE with diversity chosen.")
      }

      if (richness_type == "unique"){
        outfile_list[[i]] <- colSums(outfile_list[[i]] > 0)
      }

      if (richness_type == "cumulative"){
        outfile_list[[i]] <- colSums(t(apply(outfile_list[[i]] > 0, 1, cumsum)) > 0)
      }


      if (richness_type == "new"){
        outfile_list[[i]] <- colSums(t(apply(t(apply(outfile_list[[i]] > 0, 1, cumsum)) > 0,1, cumsum)) == 1 )
      }

      outfile_list[[i]] <- as.data.frame(outfile_list[[i]])

    }


    names(outfile_list) <- outfile_names

    merged_data <- merger(outfile_list, names = outfile_names)
    print(outfile_list)
    melty <- melt(as.matrix(merged_data))
    print(melty)
    melty <- melty[complete.cases(melty),]
    ggplot(melty, aes(x = Var1, y = value, group = Var2, colour = Var2))+
      geom_line(size = line_size)+
      geom_point(size = point_size)+
      ylab(paste0(richness_type, " Barcodes"))+
      xlab("Month")+
      coord_cartesian(ylim = c(y_lower, y_upper))+
      theme_bw()

  }

  else{
    dfr <- data.frame()
    for(i in 1:length(outfile_list)){

      outfile_list[[i]] <- outfile_list[[i]][,outfile_samples_list[[i]], drop = FALSE]
      outfile_list[[i]][outfile_list[[i]] < threshold_list[[i]]] <- 0

      listy = split(as.data.frame(t(outfile_list[[i]])), rep(1:length(celltypes_list[[i]]), length(outfile_months_list[[i]])))
      listy = lapply(listy, t)
      for(j in 1:length(listy)){
        colnames(listy[[j]]) <- outfile_months_list[[i]]
        if (richness_type == "unique"){
          listy[[j]] <- colSums(listy[[j]] > 0)
        }

        if (richness_type == "cumulative"){
          listy[[j]] <- colSums(t(apply(listy[[j]] > 0, 1, cumsum)) > 0)
        }


        if (richness_type == "new"){
          listy[[j]] <- colSums(t(apply(t(apply(listy[[j]] > 0, 1, cumsum)) > 0,1, cumsum)) == 1 )
        }

        if (richness_type == "diversity"){
          listy[[j]] <- diversity(listy[[j]], MARGIN = 2, index = diversity_index)
        }

        listy[[j]] <- as.data.frame(listy[[j]])
      }
      melty <- melt(as.matrix(merger(listy, names = paste(outfile_names[[i]], celltypes_list[[i]]))))
      melty$Monkey <- outfile_names[[i]]
      print(melty)
      dfr <- rbind(dfr, melty)

    }
    dfr <- dfr[complete.cases(dfr),]
    print(dfr)

    ggplot(dfr, aes(x = Var1, y = value, group = Var2, colour = Var2))+
      geom_line(size = line_size)+
      geom_point(size = point_size)+
      ylab(paste0(richness_type, " Barcodes"))+
      xlab("Month")+
      coord_cartesian(ylim = c(y_lower, y_upper))+
      theme_bw()




  }





}

merger <- function(listy, names = NULL, na_to_zeros = FALSE){
  df <- data.frame()
  for(i in listy){
    df <- merge(df,i, by = 0, all = TRUE)
    rownames(df) <- df$Row.names
    df$Row.names = NULL
  }

  if (na_to_zeros == TRUE){
    df[is.na(df)] <- 0
  }

  colnames(df) <- names

  return(df)
}
