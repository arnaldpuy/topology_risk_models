
# Build the static call graph, attach cyclomatic complexity (mean-filled) and
# compute path-level risk P_k with its (alpha, beta, gamma, p) ensemble --------

build_static_fun <- function(el, cc_model) {

  cc <- fread("datasets/cyclomatic_complexity_functions.csv")[model == cc_model]
  cc[, name_norm := norm_name_fun(name)]
  cc_u <- cc[, .(cyclomatic_complexity = max(cyclomatic_complexity, na.rm = TRUE)),
             by = name_norm]
  mean_cc <- mean(cc_u$cyclomatic_complexity, na.rm = TRUE)

  g   <- as_tbl_graph(el, directed = TRUE)
  nm  <- igraph::V(as.igraph(g))$name
  ccm <- cc_u[match(nm, name_norm), cyclomatic_complexity]
  ccm[is.na(ccm)] <- mean_cc
  g   <- g |> activate(nodes) |> mutate(cyclomatic_complexity = ccm)

  out <- all_paths_fun(graph = g, alpha = 0.6, beta = 0.3, gamma = 0.1, p = 1,
                       complexity_col = "cyclomatic_complexity")
  ua  <- uncertainty_fun(all_paths_out = out, N = 2^10, order = "first")

  paths <- as.data.table(out$paths)
  uap   <- as.data.table(ua$paths)
  uap[, P_k_q50 := vapply(uncertainty_analysis, median, numeric(1), na.rm = TRUE)]
  uap[, P_k_mean := vapply(uncertainty_analysis, mean, numeric(1), na.rm = TRUE)]
  uap[, P_k_sd := vapply(uncertainty_analysis, sd, numeric(1), na.rm = TRUE)]
  uap[, CV_k := P_k_sd / pmax(P_k_mean, 1e-12)]
  paths <- merge(paths, uap[, .(path_id, P_k_q50, CV_k)], by = "path_id")

  list(paths = paths, nodes = as.data.table(out$nodes))
}
