
# FUNCTION TO COMPUTE SLOPE ####################################################


slope_fun <- function(x) {
  
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) <= 1) return(0)
  
  as.numeric(coef(lm(x ~ seq_along(x)))[2])
}
