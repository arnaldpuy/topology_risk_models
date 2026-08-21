
# Join per-edge runtime call counts onto each static path: exercised flag,
# runtime intensity E_k (sum) and path-mean node metrics -----------------------

join_runtime_fun <- function(static, rt) {

  paths <- copy(static$paths)
  nodes <- static$nodes

  res <- paths[, {
    ns <- unlist(path_nodes)
    ed <- data.table(cn = norm_name_fun(ns[-length(ns)]),
                     ce = norm_name_fun(ns[-1L]))
    j  <- rt[ed, on = c("cn", "ce")]
    rc <- j$runtime_calls; ex <- !is.na(rc); rc[is.na(rc)] <- 0
    list(any_exercised = any(ex), Ek_sum = sum(rc))
  }, by = .(path_id)]
  paths <- merge(paths, res, by = "path_id")

  pm <- paths[, {
    sub <- nodes[name %in% norm_name_fun(unlist(path_nodes))]
    list(pm_cc = mean(sub$cyclomatic_complexity, na.rm = TRUE),
         pm_indeg = mean(sub$indeg, na.rm = TRUE),
         pm_btw = mean(sub$btw, na.rm = TRUE))
  }, by = .(path_id)]

  merge(paths, pm, by = "path_id")
}
