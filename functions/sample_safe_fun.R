
# FUNCTIO TO SAMPLE IN A SAFE WAY ##############################################

sample_safe_fun <- function(x, k) {
  
  x <- unique(x)
  
  if (length(x) == 0) return(x[0])
  k <- min(k, length(x))
  
  if (k <= 0) return(x[0])
  
  sample(x, k, replace = FALSE)
  
}