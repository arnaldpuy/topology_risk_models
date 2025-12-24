
# FUNCTION TO PLOT RISKY PATHS #################################################

# Define labels for better plotting --------------------------------------------

lab_expr <- c(b1 = expression(C %in% "(" * 0 * ", 10" * "]"),
              b2 = expression(C %in% "(" * 10 * ", 20" * "]"),
              b3 = expression(C %in% "(" * 20 * ", 50" * "]"),
              b4 = expression(C %in% "(" * 50 * ", " * infinity * ")"))

# Function ---------------------------------------------------------------------

plot_top_paths_fun <- function(call_g, paths_tbl, model.name = "", language = "") {
  
  if (nrow(paths_tbl) == 0) {
    message("No paths in paths_tbl; skipping plot for: ", model.name)
    return(invisible(NULL))
  }
  
  # ---- Identify top 10 risky paths -------------------------------------------
  
  top_paths <- paths_tbl %>% 
    arrange(desc(p_path_fail)) %>% 
    slice_head(n = 10)
  
  k <- min(10, nrow(top_paths))
  top_k_paths <- top_paths %>% slice_head(n = k)
  
  # ---- Build edge list -------------------------------------------------------
  
  path_edges_all <- purrr::imap_dfr(top_k_paths$path_nodes, function(nodes_vec, pid) {
    
    tibble(from = head(nodes_vec, -1), 
           to = tail(nodes_vec, -1),
           path_id = pid,
           risk_sum = top_k_paths$risk_sum[pid])
  })
  
  ig2 <- as.igraph(call_g)
  
  edge_df_names <- igraph::as_data_frame(ig2, what = "edges") %>%
    mutate(.edge_idx = dplyr::row_number())
  
  path_edges_collapsed <- path_edges_all %>%
    dplyr::group_by(from, to) %>%
    dplyr::summarise(
      path_freq      = dplyr::n(),
      risk_mean_path = mean(risk_sum, na.rm = TRUE),
      .groups = "drop"
    )
  
  edge_marks <- edge_df_names %>%
    dplyr::left_join(path_edges_collapsed, by = c("from","to")) %>%
    dplyr::mutate(on_top_path = !is.na(path_freq),
                  path_freq = ifelse(is.na(path_freq), 0L, path_freq),
                  risk_mean_path = ifelse(is.na(risk_mean_path), 0, risk_mean_path))
  
  # ---- Annotate graph --------------------------------------------------------
  
  call_g_sugi <- call_g %>%
    activate(edges) %>%
    mutate(on_top_path = edge_marks$on_top_path,
           path_freq = edge_marks$path_freq,
           risk_mean_path = edge_marks$risk_mean_path) %>%
    activate(nodes) %>%
    mutate(indeg = indeg,
           cyclo_class = case_when(cyclomatic_complexity <= 10 ~ "green",
                                   cyclomatic_complexity <= 20 ~ "orange",
                                   cyclomatic_complexity <= 50 ~ "red",
                                   cyclomatic_complexity >  50 ~ "purple",
                                   TRUE ~ "grey"), 
           complexity_category = cut(cyclomatic_complexity, 
                                     breaks = c(-Inf, 10, 20, 50, Inf),
                                     labels = c("b1","b2","b3","b4")))
  
  risky_nodes <- unique(unlist(top_k_paths$path_nodes))
  
  call_g_sugi <- call_g_sugi %>%
    activate(nodes) %>%
    mutate(on_top_node = name %in% risky_nodes)
  
  # Compute legend sizes: min, middle, max among *risky* nodes -----------------
  
  node_df <- call_g_sugi %>% 
    activate(nodes) %>% 
    as_tibble()
  
  risky_indeg <- node_df$indeg[node_df$on_top_node & is.finite(node_df$indeg)]
  risky_indeg <- sort(unique(risky_indeg))
  
  # If no risky indegree values (weird, but be safe), use all nodes ------------
  
  if (length(risky_indeg) == 0) {
    risky_indeg <- sort(unique(node_df$indeg[is.finite(node_df$indeg)]))
  }
  
  # Decide breaks depending on how many values we actually have ----------------
  
  if (length(risky_indeg) >= 3) {
    min_indeg <- min(risky_indeg)
    max_indeg <- max(risky_indeg)
    mid_indeg <- risky_indeg[ceiling(length(risky_indeg) / 2)]  # a real value from the data
    legend_breaks <- c(min_indeg, mid_indeg, max_indeg)
    
  } else {
    
    # If fewer than 3 unique values, just use what we have
    legend_breaks <- risky_indeg
  }
  
  legend_labels <- round(legend_breaks, 1)
  
  # ---- Plot ------------------------------------------------------------------
  
  p_sugi <- ggraph(call_g_sugi, layout = "sugiyama") +
    geom_edge_link0(aes(filter = !on_top_path),
                    colour = "grey80", alpha = 0.05, width = 0.3) +
    geom_edge_link0(aes(filter = on_top_path,
                        colour = risk_mean_path,
                        width = pmin(pmax(path_freq, 0.5), 3)),
                    alpha = 0.9, 
                    arrow = grid::arrow(length = unit(1, "mm"))) +
    scale_edge_colour_gradient(low = "orange", high = "red", guide = "none") +
    scale_edge_width(range = c(0.3, 2.2), guide = "none") +
    geom_node_point(size = 1.2, colour = "#BDBDBD", alpha = 0.35, show.legend = FALSE) +
    geom_node_point(aes(filter = on_top_node, size = indeg, fill = complexity_category),
                    shape = 21, alpha = 0.95, show.legend = TRUE) +
    scale_size_continuous(breaks = legend_breaks, labels = legend_labels, 
                          name = "indegree") +
    scale_fill_manual(values = c("yellowgreen", "orange", "red", "purple"),
                      labels = lab_expr,
                      name = "") +
    guides(fill = "none") +
    theme_AP() +
    labs(x = "", y = "") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "right", 
          plot.margin = margin(1, 1, 1, 1),
          panel.spacing = unit(0.1, "lines")) + 
    ggtitle(paste(model.name, ": ", language, sep = ""))
  
  print(p_sugi)
  invisible(p_sugi)
}


