extract_sa_fun <- function(dt, id_cols = c("model", "name", "risk_form")) {
  out_list <- lapply(seq_len(nrow(dt)), function(i) {
    sa <- dt$sensitivity_analysis[[i]]
    if (is.null(sa)) return(NULL)
    
    # sensobol object -> its $results tibble/data.frame
    res <- as.data.table(sa$results)
    
    # attach identifiers (model, name, risk_form, ...)
    for (col in id_cols) {
      if (col %in% names(dt)) {
        res[, (col) := dt[[col]][i] ]
      }
    }
    res
  })
  
  rbindlist(out_list, use.names = TRUE, fill = TRUE)
}