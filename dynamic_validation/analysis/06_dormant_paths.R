# 06_dormant_paths.R - Analysis 5 (dormant high-risk paths)
# Descriptive: which top-decile P_k paths are NOT exercised, and what do they
# look like? The framing argument is "static analysis identifies vulnerabilities
# runtime profiling alone misses".
#
# Top-decile is defined *within each hops bin* (mirroring Analysis 1b) so the
# dormant set is not itself a length artefact - long paths are mechanically
# less likely to have every edge exercised, and we don't want the dormancy
# claim to ride on that.
#
# Runs today on the binary edge-presence layer; no runtime counts needed.

source("dynamic_validation/analysis/00_setup.R")

RISK_FORM_PRIMARY <- "additive"
TOP_DECILE        <- 0.90
TOP_K_FUNCTIONS   <- 3L
MIN_PATHS_PER_BIN <- 20L  # mirrors 02_coverage.R

joined_files <- list.files(out_dir, pattern = "_paths_joined\\.csv$", full.names = TRUE)
if (!length(joined_files)) stop("No joined files - run 01_join_static_dynamic.R first.")

paths_df <- load_full_paths_df()
node_df  <- load_full_node_df()

dormant_summary   <- list()
dormant_paths_out <- list()

for (f in joined_files) {

  dt <- fread(f)
  m  <- unique(dt$model)
  dt <- dt[risk_form == RISK_FORM_PRIMARY & !is.na(P_k_q50)]
  if (!nrow(dt)) next

  dt[, hops_bin := cut(hops,
                       breaks = c(0, 1, 2, 4, 7, Inf),
                       labels = c("1", "2", "3-4", "5-7", "8+"),
                       right  = TRUE)]

  # Within-length top decile of P_k: high-risk *for that length class*.
  hi <- dt[
    ,
    {
      if (.N < MIN_PATHS_PER_BIN || length(unique(P_k_q50)) < 2L) {
        data.table()
      } else {
        qhi    <- quantile(P_k_q50, TOP_DECILE, na.rm = TRUE)
        cv_med <- median(CV_k, na.rm = TRUE)
        .SD[P_k_q50 >= qhi & CV_k < cv_med]
      }
    },
    by = .(model, hops_bin)
  ]

  if (!nrow(hi)) {
    cat("[dormant] ", m, ": no hops bin had enough paths for a within-length ",
        "top-decile cut. Skipping.\n", sep = "")
    next
  }

  # Looser dormancy: not even one edge exercised.
  dormant_loose  <- hi[any_exercised == FALSE | is.na(any_exercised)]
  # Stricter dormancy: not every edge exercised.
  dormant_strict <- hi[all_exercised == FALSE | is.na(all_exercised)]

  dormant_summary[[m]] <- hi[
    ,
    .(n_top_decile_in_bin = .N),
    by = .(model, hops_bin)
  ][
    merge(dormant_loose[, .N, by = .(model, hops_bin)],
          dormant_strict[, .N, by = .(model, hops_bin)],
          by = c("model", "hops_bin"),
          all = TRUE,
          suffixes = c("_loose", "_strict")),
    on = c("model", "hops_bin"),
    nomatch = NA
  ][
    , .(model, hops_bin,
        n_top_decile = n_top_decile_in_bin,
        n_dormant_loose     = fcoalesce(N_loose,  0L),
        n_dormant_strict    = fcoalesce(N_strict, 0L))
  ][
    , `:=`(
      frac_dormant_loose  = n_dormant_loose  / n_top_decile,
      frac_dormant_strict = n_dormant_strict / n_top_decile
    )
  ]

  if (!nrow(dormant_loose)) next

  per_path_nodes <- expand_paths_to_nodes(m, RISK_FORM_PRIMARY,
                                          paths_df = paths_df,
                                          node_df  = node_df)
  per_path_nodes <- per_path_nodes[path_id %in% dormant_loose$path_id]

  topk <- per_path_nodes[order(-risk_score),
                         head(.SD, TOP_K_FUNCTIONS),
                         by = .(model, path_id)]
  topk_collapsed <- topk[, .(
    top_functions = paste(name, collapse = " | "),
    top_files     = paste(unique(file), collapse = " | ")
  ), by = .(model, path_id)]

  out <- merge(dormant_loose[, .(model, path_id, hops_bin, hops, P_k_q50, CV_k)],
               topk_collapsed, by = c("model", "path_id"), all.x = TRUE)
  dormant_paths_out[[m]] <- out[order(hops_bin, -P_k_q50)]
}

if (length(dormant_summary)) {
  fwrite(rbindlist(dormant_summary, fill = TRUE),
         file.path(syn_dir, "05_dormant_summary.csv"))
  fwrite(rbindlist(dormant_paths_out, fill = TRUE),
         file.path(syn_dir, "05_dormant_paths_table.csv"))
  cat("\n[06_dormant_paths] within-length top-decile dormancy:\n")
  print(rbindlist(dormant_summary, fill = TRUE))
  cat("\n[06_dormant_paths] wrote summary + per-path table to ", syn_dir, "\n", sep = "")
} else {
  cat("\n[06_dormant_paths] no models produced dormant tables - check inputs.\n")
}
