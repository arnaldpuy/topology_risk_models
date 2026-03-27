
unnest_paths_tbl_fun <- function(dt, paths_col, slot) {
  out_list <- lapply(seq_len(nrow(dt)), function(i) {
    x <- dt[[paths_col]][[i]][[slot]] 
    if (is.null(x)) return(NULL)
    x_dt <- as.data.table(x)
    x_dt[, model := dt$model[i]]
    x_dt
  })
  
  rbindlist(out_list, use.names = TRUE, fill = TRUE)
}