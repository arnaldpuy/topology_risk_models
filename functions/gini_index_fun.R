
# FUNCTION TO COMPUTE GINI INDEX ##############################################

gini_index_fun <- function(x) {
  
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) <= 1) return(0)
  
  ineq::Gini(x)
}