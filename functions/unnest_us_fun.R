
unnest_us_fun <- function(dt, us_col, slot) {
  out_list <- lapply(seq_len(nrow(dt)), function(i) {
    x <- dt[[us_col]][[i]]      # e.g. uncertainty_sensitivity_additive[[i]]
    if (is.null(x)) return(NULL)
    
    slot_obj <- x[[slot]]       # x$nodes or x$paths
    if (is.null(slot_obj)) return(NULL)
    
    res <- as.data.table(slot_obj)
    res[, model := dt$model[i]] # keep model (add more keys if needed)
    res
  })
  
  rbindlist(out_list, use.names = TRUE, fill = TRUE)
}