custom_log <- function(x, log_choice){
  x <- log(x, log_choice)
  #sets the -inf values to be the minimum of the possible value in the data frame -1
  x[x == -Inf] <- (min(x[x > -Inf]) - 1)
  return(x)
}
