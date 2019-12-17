subset_SE <- function(your_SE, ...) {
  arguments <- list(...)
  if(length(arguments) > 0){
    for(i in 1:length(arguments)){
      subset_name <- names(arguments)[[i]]
      subset_vars <- arguments[[i]]
      your_SE <- your_SE[, your_SE[[subset_name]] %in% subset_vars]
    }
  }
  return(your_SE)
}
