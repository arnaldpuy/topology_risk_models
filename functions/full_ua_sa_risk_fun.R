
# FULL UNCERTAINTY AND SENSITIVITY ANALYSIS FUNCTION ###########################
################################################################################

# Wrapper ---------------------------------------------------------------------

risk_ua_sa_fun <- function(node_df, sample_matrix, N, params, order = order) {
  
  # Extract cyclo, in-degree and betweenness values -------------------
  
  cyclo_sc <- node_df$cyclo_sc
  indeg_sc <- node_df$indeg_sc
  btw_sc <- node_df$btw_sc
  
  # Extract only alpha, beta and gamma columns from the sample matrix---
  
  mat <- sample_matrix[, c("alpha", "beta", "gamma")]
  
  # Calculate risk score for each row ---------------------------------
  
  y <- mat[, 1] * cyclo_sc + mat[, 2] * indeg_sc + mat[, 3] * btw_sc
  
  # Calculate sobol indices -------------------------------------------
  
  ind <- sobol_indices(Y = y, N = N, params = params, order = order)
  
  # Wrap up and output ------------------------------------------------
  
  output <- list(y[1:(2*N)], ind)
  names(output) <- c("ua", "sa")
  
  return(output)
}

#################################################################################

# FULL FUNCTION ################################################################

full_ua_sa_risk_fun <- function(node_df, paths_tbl, N, order) {
  
  # Define sample matrix -------------------------------------------------------
  
  params <- c("a_raw", "b_raw", "c_raw")
  mat <- sobol_matrices(N = N, params = params, order  = order)
  s <- rowSums(mat)
  
  # Define alpha, beta, gamma --------------------------------------------------
  
  alpha <- mat[, "a_raw"] / s
  beta <- mat[, "b_raw"] / s
  gamma <- mat[, "c_raw"] / s
  
  mat <- cbind(mat, alpha, beta, gamma)
  
  # Run risk_ua_sa_fun for every row -------------------------------------------
  
  node_df <- data.table(node_df)
  out <- list()
  
  for(i in 1:nrow(node_df)) {
    
    out[[i]] <- risk_ua_sa_fun(node_df = node_df[i, ], 
                               sample_matrix = mat, 
                               N = N, 
                               params = params,
                               order = order)
    
  }
  
  # Add names of nodes ---------------------------------------------------------
  
  names(out) <- node_df$name
  
  # Retrieve the uncertainty analysis ------------------------------------------
  
  node_df[, uncertainty_analysis:= lapply(name, function(n) out[[n]][["ua"]])]
  
  # Retrieve the sensitivity analysis ------------------------------------------
  
  node_df[, sensitivity_analysis:= lapply(name, function(n) out[[n]][["sa"]])]
  node_df[, sensitivity_indices:= lapply(sensitivity_analysis, function(x) {x$results})]
  node_df[, sensitivity_analysis:= NULL]
  
  # # Define function to calculate P_k -----------------------------------------
  
  path_prob_from_nodes <- function(risks) 1 - prod(1 - risks)
  
  # Calculate P_k for every uncertainty analysis in every node -----------------
  
  P_k <- list()
  
  for (i in 1:nrow(paths_tbl)) {
    
    tmp <- node_df[name %in% paths_tbl$path_nodes[[i]]]
    ua_mat <- do.call(rbind, tmp$uncertainty_analysis)
    P_k[[i]] <- apply(ua_mat, 2, path_prob_from_nodes)
    
  }
  
  # Name with the Path ID ------------------------------------------------------
  
  names(P_k) <- paths_tbl$path_id
  
  # Retrieve data --------------------------------------------------------------
  
  paths_tbl <- data.table(paths_tbl)
  paths_tbl[, P_k_vec:= P_k[path_id]]
  
  # Compute mean, min and max to plot errorbars --------------------------------
  
  paths_tbl[, `:=`(P_k_min  = sapply(P_k_vec, min),
                   P_k_mean = sapply(P_k_vec, mean),
                   P_k_max  = sapply(P_k_vec, max), 
                   P_k_q025  = sapply(P_k_vec, quantile, probs = 0.025),
                   P_k_q50  = sapply(P_k_vec, quantile, probs = 0.5),
                   P_k_q975  = sapply(P_k_vec, quantile, probs = 0.975))]
  
  paths_tbl_sorted <- paths_tbl[order(P_k_mean)]
  
  # Return output --------------------------------------------------------------
  
  output <- list(node_df, paths_tbl_sorted)
  names(output) <- c("nodes", "paths")
  
  return(output)
}


