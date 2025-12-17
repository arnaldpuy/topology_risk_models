
# FUNCTION TO COMPUTE ALL PATHS AND CALCULATE METRICS ##########################

all_paths_fun <- function(node_df, graph) {
  
  # Identify entries and sinks -------------------------------------------------
  
  entry_ids <- which(node_df$indeg == 0)
  sink_ids  <- which(node_df$outdeg == 0)
  
  if (length(entry_ids) == 0 || length(sink_ids) == 0) {
    return(tibble())
  }
  
  # Convert to igraph ----------------------------------------------------------
  
  ig <- as.igraph(graph)
  
  # Distance matrix ------------------------------------------------------------
  
  dist_mat <- igraph::distances(ig, v = entry_ids, to = sink_ids, mode = "out")
  
  # Build all entry–sink index pairs -------------------------------------------
  
  pairs_st <- expand.grid(s = seq_along(entry_ids), t = seq_along(sink_ids))
  
  # Exclude cases where entry and sink are the *same vertex* -------------------
  
  pairs_st <- subset(pairs_st, entry_ids[s] != sink_ids[t])
  
  # Keep only reachable pairs (finite distance) --------------------------------
  
  pairs_st <- subset(pairs_st, dist_mat[cbind(s, t)] < Inf)
  
  if (nrow(pairs_st) == 0) {
    return(tibble())
  }
  
  # Compute all simple paths (nested list) -------------------------------------
  
  all_paths_nested <- lapply(seq_len(nrow(pairs_st)), function(i) {
    from_v <- entry_ids[pairs_st$s[i]]
    to_v  <- sink_ids[pairs_st$t[i]]
    igraph::all_simple_paths(ig, from = from_v, to = to_v, mode = "out")
  })
  
  # Flatten: list-of-lists to list of paths ------------------------------------
  
  all_paths <- purrr::flatten(all_paths_nested)
  if (length(all_paths) == 0) return(tibble())
  
  # Create path data frame -----------------------------------------------------
  
  paths_tbl <- map_dfr(seq_along(all_paths), function(k) {
    v <- all_paths[[k]]
    v_names <- names(v)
    ps <- node_df$risk_score[match(v_names, node_df$name)]
    
    # Probability that at least one node fails on the path 
    # (union of independent events) --------------------------------------------
    
    Pk <- if (length(ps) == 0) 0 else 1 - prod(1 - ps)
    
    tibble(path_id = k,
           path_nodes = list(v_names),
           path_str = paste(v_names, collapse = " \u2192 "),
           hops = length(v_names) - 1,
           p_path_fail = Pk,
           gini_node_risk = gini_index_fun(ps),
           risk_slope = slope_fun(ps),
           risk_mean = mean(ps, na.rm = TRUE), 
           risk_sum = sum(ps, na.rm = TRUE))
  })
  
  return(paths_tbl)
}

