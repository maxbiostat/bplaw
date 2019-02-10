compress_data <- function(x){
  raw <- table(x)
  vv <- as.numeric(names(raw))
  ns <- as.numeric(raw)
  return(list(v = vv, fs = ns, K = length(vv)))
}