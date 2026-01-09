

sensobol::load_packages(c("data.table", "tidyverse", "openxlsx", "scales", 
                          "cowplot", "readxl", "ggrepel", "tidytext", "here", 
                          "tidygraph", "igraph", "foreach", "parallel", "ggraph", 
                          "tools", "purrr", "sensobol", "benchmarkme"))

# Create custom theme ----------------------------------------------------------

theme_AP <- function() {
  theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.key = element_rect(fill = "transparent", color = NA), 
          strip.background = element_rect(fill = "white"), 
          legend.text = element_text(size = 7.3), 
          axis.title = element_text(size = 10),
          legend.key.width = unit(0.4, "cm"), 
          legend.key.height = unit(0.4, "cm"), 
          legend.key.spacing.y = unit(0, "lines"),
          legend.box.spacing = unit(0, "pt"),
          legend.title = element_text(size = 7.3), 
          axis.text.x = element_text(size = 7), 
          axis.text.y = element_text(size = 7), 
          axis.title.x = element_text(size = 7.3), 
          axis.title.y = element_text(size = 7.3),
          plot.title = element_text(size = 8),
          strip.text.x = element_text(size = 7.4), 
          strip.text.y = element_text(size = 7.4)) 
}


# ============================================================
# Synthetic scalability experiment (layered DAG only)
# - Sample execution paths via random walks (no exhaustive enumeration)
# - Path risk via saturating aggregation: 1 - prod(1 - r_i)
# - UA/SA via simplex sampling of (alpha,beta,gamma) and top-K stability
# - Parallelized across configurations with foreach/doParallel
# ============================================================

# ----------------------------
# 1) Synthetic call-graph generator (layered DAG)
# ----------------------------
make_synthetic_callgraph <- function(depth = 20, width = 200, branching = 3, seed = 1) {
  set.seed(seed)
  
  entry_id <- 1L
  layers <- vector("list", depth + 1L)
  layers[[1L]] <- entry_id
  
  id <- entry_id
  for (d in 2:(depth + 1L)) {
    layers[[d]] <- (id + 1L):(id + width)
    id <- id + width
  }
  exit_ids <- layers[[depth + 1L]]
  
  edges <- vector("list", depth)
  for (d in 1:depth) {
    from <- layers[[d]]
    to   <- layers[[d + 1L]]
    
    k <- pmax(1L, rpois(length(from), lambda = branching))
    ed <- lapply(seq_along(from), function(i) {
      cbind(from[i], sample(to, size = k[i], replace = TRUE))
    })
    edges[[d]] <- do.call(rbind, ed)
  }
  
  edge_mat <- do.call(rbind, edges)
  
  # explicit vertex names as characters
  edge_mat_chr <- matrix(as.character(edge_mat), ncol = 2)
  g <- graph_from_edgelist(edge_mat_chr, directed = TRUE)
  
  # clean up
  g <- delete_edges(g, E(g)[which_loop(g)])
  g <- delete_edges(g, E(g)[which_multiple(g)])
  
  # ensure all vertices exist
  all_ids <- as.character(seq_len(id))
  missing <- setdiff(all_ids, V(g)$name)
  if (length(missing) > 0) {
    g <- add_vertices(g, nv = length(missing), name = missing)
  }
  
  list(
    g = g,
    entry_name = "1",
    exit_names = as.character(exit_ids)
  )
}

# ----------------------------
# 2) Node attributes: cyclo + indeg + betweenness; rescale() normalization
#    Robust betweenness across igraph versions: estimate -> cutoff -> exact
# ----------------------------
assign_node_attributes <- function(g, seed = 1,
                                   cyclo_dist = c("lognormal", "gamma"),
                                   cyclo_scale = 1.0,
                                   btw_method = c("estimate", "cutoff", "exact"),
                                   btw_cutoff = 6) {
  set.seed(seed)
  cyclo_dist <- match.arg(cyclo_dist)
  btw_method <- match.arg(btw_method)
  
  n <- vcount(g)
  
  cyclo_raw <- switch(
    cyclo_dist,
    lognormal = rlnorm(n, meanlog = log(3), sdlog = 0.9),
    gamma     = rgamma(n, shape = 2.0, rate = 0.5)
  )
  cyclo_raw <- pmax(1, cyclo_raw * cyclo_scale)
  
  indeg_raw <- degree(g, mode = "in")
  
  btw_raw <- NULL
  if (btw_method == "estimate" &&
      exists("estimate_betweenness", where = asNamespace("igraph"), inherits = FALSE)) {
    btw_raw <- tryCatch(
      igraph::estimate_betweenness(g, directed = TRUE),
      error = function(e) NULL
    )
  }
  if (is.null(btw_raw) && btw_method %in% c("estimate", "cutoff")) {
    btw_raw <- igraph::betweenness(g, directed = TRUE, normalized = FALSE, cutoff = btw_cutoff)
  }
  if (is.null(btw_raw)) {
    btw_raw <- igraph::betweenness(g, directed = TRUE, normalized = FALSE)
  }
  
  V(g)$cyclo_raw <- as.numeric(cyclo_raw)
  V(g)$indeg_raw <- as.numeric(indeg_raw)
  V(g)$btw_raw   <- as.numeric(btw_raw)
  
  V(g)$cyclo <- as.numeric(scales::rescale(V(g)$cyclo_raw, to = c(0, 1)))
  V(g)$indeg <- as.numeric(scales::rescale(V(g)$indeg_raw, to = c(0, 1)))
  V(g)$btw   <- as.numeric(scales::rescale(V(g)$btw_raw,   to = c(0, 1)))
  
  g
}

# ----------------------------
# 3) Fast path sampling (random walks) using adjacency lists
# ----------------------------
sample_paths_fast <- function(g, entry_name, exit_names,
                              n_paths = 5000, max_steps = 500, seed = 1) {
  set.seed(seed)
  
  entry_vid <- which(V(g)$name == entry_name)
  exit_vids <- which(V(g)$name %in% exit_names)
  
  if (length(entry_vid) != 1) stop("Entry node not found uniquely by name.")
  if (length(exit_vids) == 0) stop("No exit nodes found by name.")
  
  adj <- lapply(as_adj_list(g, mode = "out"), as.integer)
  is_exit <- rep(FALSE, vcount(g))
  is_exit[exit_vids] <- TRUE
  
  paths <- vector("list", n_paths)
  ok <- logical(n_paths)
  
  for (p in seq_len(n_paths)) {
    cur <- entry_vid
    path <- integer(1 + max_steps)
    path[1L] <- cur
    len <- 1L
    steps <- 0L
    
    while (!is_exit[cur] && steps < max_steps) {
      nbrs <- adj[[cur]]
      if (length(nbrs) == 0L) break
      cur <- nbrs[sample.int(length(nbrs), 1L)]
      steps <- steps + 1L
      len <- len + 1L
      path[len] <- cur
    }
    
    if (is_exit[cur]) {
      ok[p] <- TRUE
      paths[[p]] <- path[seq_len(len)]
    } else {
      paths[[p]] <- NULL
    }
  }
  
  paths[ok]
}

# ----------------------------
# 4) Path risk scoring: saturating aggregation
#    P_k = 1 - prod(1 - r_i)
# ----------------------------
score_paths_saturating <- function(g, paths, alpha = 1/3, beta = 1/3, gamma = 1/3) {
  r_node <- alpha * V(g)$cyclo + beta * V(g)$indeg + gamma * V(g)$btw
  vapply(paths, function(p) 1 - prod(1 - r_node[p]), numeric(1))
}

# ----------------------------
# 5) UA/SA: sample (alpha,beta,gamma) on simplex; compute top-K stability
# ----------------------------
ua_sa_path_stability <- function(g, paths, K = 10, n_weight_samples = 200,
                                 weight_sampler = c("dirichlet", "uniform"),
                                 seed = 1) {
  set.seed(seed)
  weight_sampler <- match.arg(weight_sampler)
  
  W <- switch(weight_sampler,
              dirichlet = {
                X <- matrix(rexp(n_weight_samples * 3), ncol = 3)
                X / rowSums(X)
              },
              uniform = {
                X <- matrix(runif(n_weight_samples * 3), ncol = 3)
                X / rowSums(X)
              }
  )
  colnames(W) <- c("alpha", "beta", "gamma")
  
  if (length(paths) == 0) {
    return(list(weights = W, freq_topK = numeric(0), stable_paths = integer(0)))
  }
  
  topK <- vector("list", n_weight_samples)
  for (i in seq_len(n_weight_samples)) {
    sc <- score_paths_saturating(g, paths, alpha = W[i,1], beta = W[i,2], gamma = W[i,3])
    topK[[i]] <- order(sc, decreasing = TRUE)[seq_len(min(K, length(sc)))]
  }
  
  counts <- integer(length(paths))
  for (i in seq_len(n_weight_samples)) counts[topK[[i]]] <- counts[topK[[i]]] + 1L
  freq <- counts / n_weight_samples
  
  stable_rank <- order(freq, decreasing = TRUE)
  
  list(
    weights = W,
    freq_topK = freq,
    stable_paths = stable_rank[seq_len(min(K, length(stable_rank)))]
  )
}

# ----------------------------
# 6) Parallel scalability runner (foreach)
# ----------------------------
run_scalability_experiment_parallel <- function(configs,
                                                n_paths = 5000,
                                                n_weight_samples = 200,
                                                K = 10,
                                                seed = 42,
                                                n_cores = max(1L, parallel::detectCores() - 1L),
                                                btw_method = c("estimate", "cutoff", "exact"),
                                                btw_cutoff = 6) {
  btw_method <- match.arg(btw_method)
  
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  res <- foreach(
    i = seq_along(configs),
    .combine = rbind,
    .packages = c("igraph", "scales"),
    .export = c(
      "make_synthetic_callgraph",
      "assign_node_attributes",
      "sample_paths_fast",
      "score_paths_saturating",
      "ua_sa_path_stability"
    )
  ) %dopar% {
    
    cfg <- configs[[i]]
    
    # ---- Build graph + attributes ----
    t0 <- proc.time()[["elapsed"]]
    cg <- make_synthetic_callgraph(
      depth = cfg$depth, width = cfg$width, branching = cfg$branching, seed = seed + i
    )
    g <- assign_node_attributes(
      cg$g,
      seed = seed + 1000 + i,
      btw_method = btw_method,
      btw_cutoff = btw_cutoff
    )
    t_build <- proc.time()[["elapsed"]] - t0
    
    # ---- Sample paths ----
    t1 <- proc.time()[["elapsed"]]
    paths <- sample_paths_fast(
      g, cg$entry_name, cg$exit_names,
      n_paths = n_paths,
      max_steps = max(500, cfg$depth * 5),
      seed = seed + 2000 + i
    )
    t_paths <- proc.time()[["elapsed"]] - t1
    
    # ---- UA/SA ----
    t2 <- proc.time()[["elapsed"]]
    uasa <- ua_sa_path_stability(
      g, paths,
      K = K,
      n_weight_samples = n_weight_samples,
      seed = seed + 3000 + i
    )
    t_uasa <- proc.time()[["elapsed"]] - t2
    
    # Stability summaries
    if (length(uasa$freq_topK) == 0) {
      top1_stability <- NA_real_
      topK_median_stability <- NA_real_
      topK_mean_stability <- NA_real_
    } else {
      stable_freqs <- sort(uasa$freq_topK, decreasing = TRUE)
      top1_stability <- stable_freqs[1L]
      topK_vec <- stable_freqs[seq_len(min(K, length(stable_freqs)))]
      topK_median_stability <- median(topK_vec)
      topK_mean_stability <- mean(topK_vec)
    }
    
    data.frame(
      depth = cfg$depth,
      width = cfg$width,
      branching = cfg$branching,
      n_nodes = vcount(g),
      n_edges = ecount(g),
      n_sampled_paths = length(paths),
      time_build = t_build,
      time_sample_paths = t_paths,
      time_uasa = t_uasa,
      top1_stability = top1_stability,
      topK_median_stability = topK_median_stability,
      topK_mean_stability = topK_mean_stability,
      btw_method = btw_method,
      btw_cutoff = btw_cutoff,
      stringsAsFactors = FALSE
    )
  }
  
  res
}

# ----------------------------
# Example configs (intermediate sizes)
# ----------------------------
configs <- list(
  list(depth = 20, width = 50,  branching = 2.5),
  list(depth = 20, width = 100, branching = 2.5),
  list(depth = 20, width = 150, branching = 2.5),
  list(depth = 20, width = 200, branching = 2.5),
  list(depth = 40, width = 200, branching = 3.0),
  list(depth = 50, width = 250, branching = 3.0),
  list(depth = 60, width = 300, branching = 3.0)
)

# ----------------------------
# Run
# ----------------------------
res <- run_scalability_experiment_parallel(
  configs = configs,
  n_paths = 2000,
  n_weight_samples = 2000,
  K = 10,
  seed = 42,
  n_cores = 8,
  btw_method = "estimate",  # falls back to cutoff/exact if unavailable
  btw_cutoff = 6
)

print(res)


###############
###############

library(ggplot2)
library(tidyr)
library(dplyr)


# Time vs n_nodes --------------------------------------------------------------

res_time_long <- res %>%
  select(n_nodes, time_build, time_sample_paths, time_uasa) %>%
  pivot_longer(cols = starts_with("time_"),
               names_to = "stage",
               values_to = "seconds") %>%
  mutate(stage = recode(stage,
                        time_build = "Graph + attributes",
                        time_sample_paths = "Path sampling",
                        time_uasa = "Uncertainty analysis"))

p_time <- ggplot(res_time_long, aes(x = n_nodes, y = seconds, group = stage, color = stage)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = "Nº nodes", y = "Time (s)") +
  scale_color_discrete(name = "") +
  theme_AP() +
  theme(legend.position = c(0.7,0.4))

p_time


# Stability vs n_nodes ---------------------------------------------------------

res_stab_long <- res %>%
  select(n_nodes, top1_stability, topK_median_stability, topK_mean_stability) %>%
  pivot_longer(cols = c(top1_stability, topK_median_stability, topK_mean_stability),
               names_to = "metric", values_to = "stability") %>%
  mutate(metric = recode(metric,
                         top1_stability = "Top-1 stability",
                         topK_median_stability = "Top-10 median stability",
                         topK_mean_stability = "Top-10 mean stability"))

p_stab <- ggplot(res_stab_long, aes(x = n_nodes, y = stability, group = metric, color = metric)) +
  geom_line() +
  geom_point() +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Nº nodes", y = "Stability", shape = "Metric") +
  theme_AP() +
  scale_color_manual(values = c("blue", "darkgreen", "orange"), 
                     name = "") +
  theme(legend.position = c(0.4, 0.22))

p_stab


# Top-K median stability vs n_edges---------------------------------------------

p_stab_edges <- ggplot(res, aes(x = n_edges, y = topK_median_stability)) +
  geom_line() +
  geom_point() +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(breaks = breaks_pretty(n = 3)) +
  labs(x = "Nº edges", y = "Top-10 med. stability") +
  theme_AP()

p_stab_edges


# Stability vs depth (architecture) --------------------------------------------

p_stab_depth <- ggplot(res, aes(x = depth, y = top1_stability)) +
  geom_line() +
  geom_point() +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Depth", y = "Top-1 stability") +
  theme_AP()

p_stab_depth

# Merge ########################################################################

top <- plot_grid(p_time, p_stab, ncol = 2, labels = "auto")
bottom <- plot_grid(p_stab_depth, p_stab_edges, ncol = 2, labels = c("c", "d"))
plot_grid(top, bottom, ncol = 1)


################################################################################
################################################################################


suppressPackageStartupMessages({
  library(igraph)
  library(scales)
  library(foreach)
  library(doParallel)
})

make_synthetic_callgraph <- function(depth = 20, width = 200, branching = 3, seed = 1) {
  set.seed(seed)
  
  entry_id <- 1L
  layers <- vector("list", depth + 1L)
  layers[[1L]] <- entry_id
  
  id <- entry_id
  for (d in 2:(depth + 1L)) {
    layers[[d]] <- (id + 1L):(id + width)
    id <- id + width
  }
  exit_ids <- layers[[depth + 1L]]
  
  edges <- vector("list", depth)
  for (d in 1:depth) {
    from <- layers[[d]]
    to   <- layers[[d + 1L]]
    
    k <- pmax(1L, rpois(length(from), lambda = branching))
    ed <- lapply(seq_along(from), function(i) {
      cbind(from[i], sample(to, size = k[i], replace = TRUE))
    })
    edges[[d]] <- do.call(rbind, ed)
  }
  
  edge_mat <- do.call(rbind, edges)
  
  # IMPORTANT FIX: make vertex names explicit
  edge_mat_chr <- matrix(as.character(edge_mat), ncol = 2)
  
  g <- graph_from_edgelist(edge_mat_chr, directed = TRUE)
  
  # clean up
  g <- delete_edges(g, E(g)[which_loop(g)])
  g <- delete_edges(g, E(g)[which_multiple(g)])
  
  # If you *really* want to ensure all vertices exist, do it safely:
  all_ids <- as.character(seq_len(id))
  missing <- setdiff(all_ids, V(g)$name)
  if (length(missing) > 0) {
    g <- add_vertices(g, nv = length(missing), name = missing)
  }
  
  list(
    g = g,
    entry_name = "1",
    exit_names = as.character(exit_ids)
  )
}


# ----------------------------
# 2) Node attributes: cyclo (right tail) + indeg + btw; rescale() normalization
# ----------------------------
assign_node_attributes <- function(g, seed = 1,
                                   cyclo_dist = c("lognormal", "gamma"),
                                   cyclo_scale = 1.0) {
  set.seed(seed)
  cyclo_dist <- match.arg(cyclo_dist)
  
  n <- vcount(g)
  
  cyclo_raw <- switch(
    cyclo_dist,
    lognormal = rlnorm(n, meanlog = log(3), sdlog = 0.9),
    gamma     = rgamma(n, shape = 2.0, rate = 0.5)
  )
  cyclo_raw <- pmax(1, cyclo_raw * cyclo_scale)
  
  indeg_raw <- degree(g, mode = "in")
  btw_raw   <- betweenness(g, directed = TRUE, normalized = FALSE)
  
  V(g)$cyclo_raw <- as.numeric(cyclo_raw)
  V(g)$indeg_raw <- as.numeric(indeg_raw)
  V(g)$btw_raw   <- as.numeric(btw_raw)
  
  # Rescale to [0,1] (risk computed on rescaled values)
  V(g)$cyclo <- as.numeric(scales::rescale(V(g)$cyclo_raw, to = c(0, 1)))
  V(g)$indeg <- as.numeric(scales::rescale(V(g)$indeg_raw, to = c(0, 1)))
  V(g)$btw   <- as.numeric(scales::rescale(V(g)$btw_raw,   to = c(0, 1)))
  
  g
}

# ----------------------------
# 3) Sample execution paths (random walks from entry to any exit)
# IMPORTANT: convert entry/exit NAMES -> vertex indices.
# ----------------------------
sample_paths <- function(g, entry_name, exit_names, n_paths = 5000, max_steps = 500, seed = 1) {
  set.seed(seed)
  
  entry_vid <- which(V(g)$name == entry_name)
  exit_vids <- which(V(g)$name %in% exit_names)
  
  if (length(entry_vid) != 1) stop("Entry node not found uniquely by name.")
  if (length(exit_vids) == 0) stop("No exit nodes found by name.")
  
  paths <- vector("list", n_paths)
  
  for (p in seq_len(n_paths)) {
    cur <- entry_vid
    path <- cur
    steps <- 0L
    
    while (!(cur %in% exit_vids) && steps < max_steps) {
      nbrs <- neighbors(g, cur, mode = "out")
      if (length(nbrs) == 0L) break
      cur <- as.integer(sample(nbrs, 1))
      path <- c(path, cur)
      steps <- steps + 1L
    }
    paths[[p]] <- path
  }
  
  # Keep only those that reached an exit
  ok <- vapply(paths, function(x) tail(x, 1) %in% exit_vids, logical(1))
  paths[ok]
}

# ----------------------------
# 4) Risk scoring: rescaled node risk + saturating path aggregator
# P_k = 1 - prod(1 - r_i)
# ----------------------------
score_paths_saturating <- function(g, paths, alpha = 1/3, beta = 1/3, gamma = 1/3) {
  r_node <- alpha * V(g)$cyclo + beta * V(g)$indeg + gamma * V(g)$btw
  vapply(paths, function(p) 1 - prod(1 - r_node[p]), numeric(1))
}

# ----------------------------
# 5) UA/SA: sample (alpha,beta,gamma) on simplex and compute top-K stability
# ----------------------------
ua_sa_path_stability <- function(g, paths, K = 10, n_weight_samples = 200,
                                 weight_sampler = c("dirichlet", "uniform"),
                                 seed = 1) {
  set.seed(seed)
  weight_sampler <- match.arg(weight_sampler)
  
  W <- switch(weight_sampler,
              dirichlet = {
                X <- matrix(rexp(n_weight_samples * 3), ncol = 3)
                X / rowSums(X)
              },
              uniform = {
                X <- matrix(runif(n_weight_samples * 3), ncol = 3)
                X / rowSums(X)
              }
  )
  colnames(W) <- c("alpha", "beta", "gamma")
  
  if (length(paths) == 0) {
    return(list(weights = W, freq_topK = numeric(0), stable_paths = integer(0)))
  }
  
  topK <- vector("list", n_weight_samples)
  for (i in seq_len(n_weight_samples)) {
    sc <- score_paths_saturating(g, paths, alpha = W[i,1], beta = W[i,2], gamma = W[i,3])
    topK[[i]] <- order(sc, decreasing = TRUE)[seq_len(min(K, length(sc)))]
  }
  
  counts <- integer(length(paths))
  for (i in seq_len(n_weight_samples)) counts[topK[[i]]] <- counts[topK[[i]]] + 1L
  freq <- counts / n_weight_samples
  
  stable_rank <- order(freq, decreasing = TRUE)
  
  list(
    weights = W,
    freq_topK = freq,
    stable_paths = stable_rank[seq_len(min(K, length(stable_rank)))]
  )
}

# ----------------------------
# 6) Parallel scalability runner (foreach)
# ----------------------------
run_scalability_experiment_parallel <- function(configs,
                                                n_paths = 5000,
                                                n_weight_samples = 200,
                                                K = 10,
                                                seed = 42,
                                                n_cores = max(1L, parallel::detectCores() - 1L)) {
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  res <- foreach(
    i = seq_along(configs),
    .combine = rbind,
    .packages = c("igraph", "scales"),
    .export = c(
      "make_synthetic_callgraph",
      "assign_node_attributes",
      "sample_paths",
      "score_paths_saturating",
      "ua_sa_path_stability"
    )
  ) %dopar% {
                   
                   cfg <- configs[[i]]
                   
                   # Build graph + attributes
                   t0 <- proc.time()[["elapsed"]]
                   cg <- make_synthetic_callgraph(depth = cfg$depth, width = cfg$width,
                                                  branching = cfg$branching, seed = seed + i)
                   g <- assign_node_attributes(cg$g, seed = seed + 1000 + i)
                   t_build <- proc.time()[["elapsed"]] - t0
                   
                   # Sample paths
                   t1 <- proc.time()[["elapsed"]]
                   paths <- sample_paths(g, cg$entry_name, cg$exit_names,
                                         n_paths = n_paths, max_steps = max(500, cfg$depth * 5),
                                         seed = seed + 2000 + i)
                   t_paths <- proc.time()[["elapsed"]] - t1
                   
                   # UA/SA stability
                   t2 <- proc.time()[["elapsed"]]
                   uasa <- ua_sa_path_stability(g, paths, K = K, n_weight_samples = n_weight_samples,
                                                seed = seed + 3000 + i)
                   t_uasa <- proc.time()[["elapsed"]] - t2
                   
                   # Summarize stability (median frequency among the most stable K paths)
                   if (length(uasa$freq_topK) == 0) {
                     topK_median_stability <- NA_real_
                   } else {
                     stable_freqs <- sort(uasa$freq_topK, decreasing = TRUE)
                     topK_median_stability <- median(stable_freqs[seq_len(min(K, length(stable_freqs)))])
                   }
                   
                   data.frame(
                     depth = cfg$depth,
                     width = cfg$width,
                     branching = cfg$branching,
                     n_nodes = vcount(g),
                     n_edges = ecount(g),
                     n_sampled_paths = length(paths),
                     time_build = t_build,
                     time_sample_paths = t_paths,
                     time_uasa = t_uasa,
                     topK_median_stability = topK_median_stability
                   )
                 }
  
  res
}

# ----------------------------
# Example configs + run
# ----------------------------
configs <- list(
  list(depth = 20, width = 50,  branching = 2.5),
  list(depth = 20, width = 200, branching = 2.5),
  list(depth = 40, width = 200, branching = 3.0),
  list(depth = 60, width = 300, branching = 3.0)
)

res <- run_scalability_experiment_parallel(
  configs = configs,
  n_paths = 5000,
  n_weight_samples = 200,
  K = 10,
  seed = 42,
  n_cores = 8
)

print(res)



# Optional: quick plots (base R; no fixed colors)
# Runtime scaling
plot(res$n_nodes, res$time_uasa, log = "xy",
     xlab = "Number of nodes (log)", ylab = "UA/SA time (s, log)",
     main = "Synthetic scalability: UA/SA runtime")
# Stability indicator
plot(res$n_nodes, res$topK_median_stability,
     xlab = "Number of nodes", ylab = "Median top-K stability",
     main = "Top-risk path stability under weight uncertainty")















# Synthetic call-graph generator for softwareRisk scalability tests
# ---------------------------------------------------------------
# Produces directed acyclic call graphs (DAGs) as tidygraph::tbl_graph
# with a node attribute `cyclomatic_complexity` (default) as expected
# by softwareRisk::all_paths_fun(..., complexity_col = "cyclomatic_complexity").
#
# You can generate a *suite* of n graphs with increasing depth/width/branching.
#
# Dependencies: tidygraph, igraph, dplyr

suppressPackageStartupMessages({
  library(igraph)
  library(tidygraph)
  library(dplyr)
})

# ---- Safe loop/multi-edge removal (igraph-version agnostic) ----
.remove_loops_and_multiedges <- function(g) {
  if (ecount(g) == 0) return(g)
  g <- delete_edges(g, E(g)[which_loop(g)])
  g <- delete_edges(g, E(g)[which_multiple(g)])
  g
}

# ---- 1) Layered DAG generator (best for controlled path growth) ----
# depth: number of transitions (layers) from sources to sinks
# width: nodes per intermediate layer (except sources/sinks if you choose)
# n_sources, n_sinks: number of entry and exit nodes
# branching: expected out-degree per node (Poisson around this, min 1 where applicable)
# p_skip: probability of adding a "skip edge" from layer l to l+2 (kept acyclic)
# p_drop: probability to drop an otherwise created edge (sparseness control)
# complexity_dist: how to sample cyclomatic complexity per node ("lognormal" or "gamma")
# complexity_scale: multiplies complexity (useful to increase mean with size)
make_synthetic_call_dag <- function(
    depth = 20,
    width = 200,
    n_sources = 1,
    n_sinks = 50,
    branching = 3,
    p_skip = 0.05,
    p_drop = 0.0,
    complexity_dist = c("lognormal", "gamma"),
    complexity_scale = 1.0,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  complexity_dist <- match.arg(complexity_dist)
  
  stopifnot(depth >= 2, width >= 1, n_sources >= 1, n_sinks >= 1, branching > 0)
  
  # Build layers: L0 sources, L1..L(depth-1) intermediate, Ldepth sinks
  layers <- vector("list", depth + 1L)
  id <- 0L
  layers[[1L]] <- (id + 1L):(id + n_sources); id <- id + n_sources
  
  for (l in 2:depth) {
    layers[[l]] <- (id + 1L):(id + width); id <- id + width
  }
  layers[[depth + 1L]] <- (id + 1L):(id + n_sinks); id <- id + n_sinks
  
  # Edges between consecutive layers
  edges <- list()
  
  for (l in 1:depth) {
    from <- layers[[l]]
    to   <- layers[[l + 1L]]
    
    # Outgoing edges per node (min 1 except we still allow dropping later)
    k <- pmax(1L, rpois(length(from), lambda = branching))
    
    el <- lapply(seq_along(from), function(i) {
      cbind(from[i], sample(to, size = k[i], replace = TRUE))
    })
    edges[[length(edges) + 1L]] <- do.call(rbind, el)
    
    # Optional skip edges to increase "architectural shortcuts" but remain acyclic
    if (p_skip > 0 && l <= depth - 1L) {
      to2 <- layers[[l + 2L]]
      # Each "from" node gets a skip edge with probability p_skip
      add_skip <- runif(length(from)) < p_skip
      if (any(add_skip)) {
        skip_mat <- cbind(from[add_skip], sample(to2, size = sum(add_skip), replace = TRUE))
        edges[[length(edges) + 1L]] <- skip_mat
      }
    }
  }
  
  edge_mat <- do.call(rbind, edges)
  
  # Optional random edge dropping to control density/path explosion
  if (p_drop > 0) {
    keep <- runif(nrow(edge_mat)) > p_drop
    edge_mat <- edge_mat[keep, , drop = FALSE]
    if (nrow(edge_mat) == 0) {
      stop("All edges dropped; reduce p_drop.")
    }
  }
  
  g <- graph_from_edgelist(edge_mat, directed = TRUE)
  g <- .remove_loops_and_multiedges(g)
  
  # Ensure all nodes are present (igraph will drop isolates otherwise)
  # Add any missing vertices explicitly
  missing <- setdiff(seq_len(id), as.integer(V(g)$name))
  if (length(missing) > 0) {
    g <- add_vertices(g, length(missing), name = as.character(missing))
  }
  
  # Convert vertex names to integers for convenience
  V(g)$id_int <- as.integer(V(g)$name)
  
  # Ensure required igraph vertex attribute `name` exists (softwareRisk requirement)
  # (Keep as character labels)
  V(g)$name <- as.character(V(g)$id_int)
  
  # Assign synthetic cyclomatic complexity (positive integer-ish)
  n <- vcount(g)
  cc <- switch(
    complexity_dist,
    lognormal = rlnorm(n, meanlog = log(3), sdlog = 0.6),
    gamma     = rgamma(n, shape = 2.5, rate = 0.6)
  )
  cc <- pmax(1, round(cc * complexity_scale))
  V(g)$cyclomatic_complexity <- cc
  
  # Return as tidygraph::tbl_graph; keep `name` attribute intact
  tg <- as_tbl_graph(g) %>%
    activate(nodes) %>%
    mutate(
      cyclomatic_complexity = .data$cyclomatic_complexity
    )
  
  tg
}

# ---- 2) Generate a suite of n graphs with increasing complexity ----
# Default schedule ramps depth, width, branching, and (optionally) complexity_scale.
generate_scalability_suite <- function(
    n = 10,
    seed = 1,
    depth_seq = round(seq(10, 60, length.out = n)),
    width_seq = round(seq(50, 400, length.out = n)),
    branching_seq = seq(2.0, 4.0, length.out = n),
    n_sources = 1,
    n_sinks_seq = round(seq(20, 200, length.out = n)),
    p_skip = 0.05,
    p_drop = 0.0,
    complexity_scale_seq = seq(1.0, 1.0, length.out = n),
    complexity_dist = "lognormal"
) {
  stopifnot(length(depth_seq) == n,
            length(width_seq) == n,
            length(branching_seq) == n,
            length(n_sinks_seq) == n,
            length(complexity_scale_seq) == n)
  
  suite <- vector("list", n)
  
  for (i in seq_len(n)) {
    tg <- make_synthetic_call_dag(
      depth = depth_seq[i],
      width = width_seq[i],
      n_sources = n_sources,
      n_sinks = n_sinks_seq[i],
      branching = branching_seq[i],
      p_skip = p_skip,
      p_drop = p_drop,
      complexity_dist = complexity_dist,
      complexity_scale = complexity_scale_seq[i],
      seed = seed + i
    )
    
    suite[[i]] <- list(
      sim_id = i,
      graph = tg,
      params = list(
        depth = depth_seq[i],
        width = width_seq[i],
        branching = branching_seq[i],
        n_sources = n_sources,
        n_sinks = n_sinks_seq[i],
        p_skip = p_skip,
        p_drop = p_drop,
        complexity_dist = complexity_dist,
        complexity_scale = complexity_scale_seq[i]
      ),
      n_nodes = igraph::vcount(tg),
      n_edges = igraph::ecount(tg)
    )
  }
  
  suite
}


# ---- 3) Convenience: suite summary table ----
suite_summary <- function(suite) {
  do.call(rbind, lapply(suite, function(s) {
    data.frame(
      sim_id = s$sim_id,
      depth = s$params$depth,
      width = s$params$width,
      branching = s$params$branching,
      n_sources = s$params$n_sources,
      n_sinks = s$params$n_sinks,
      p_skip = s$params$p_skip,
      p_drop = s$params$p_drop,
      n_nodes = s$n_nodes,
      n_edges = s$n_edges
    )
  }))
}

# -----------------------
# Example usage
# -----------------------
suite <- generate_scalability_suite(n = 10, seed = 42)
# print(suite_summary(suite))
#
# Then run your functions, e.g.:
for (s in suite) {
   g <- s$graph
   out <- softwareRisk::all_paths_fun(g, alpha = 1/3, beta = 1/3, gamma = 1/3,
                                   complexity_col = "cyclomatic_complexity")
}

