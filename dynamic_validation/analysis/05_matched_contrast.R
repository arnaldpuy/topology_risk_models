# 05_matched_contrast.R - Analysis 4 (matched-frequency contrast)
# Cleanest test of "static topology adds information beyond execution frequency".
# For each high-P_k path that is exercised, find a low-P_k path with comparable
# E_k, path length, and total statement count. Then compare static metrics on
# matched pairs.
#
# BLOCKED ON RUNTIME COUNTS: matching covariate E_k needs runtime_calls. The
# script will skip per model until counts are available.
#
# Caveat: with small samples, matching may not yield enough pairs. We mark
# results as exploratory if N pairs < 15.

source("dynamic_validation/analysis/00_setup.R")

RISK_FORM_PRIMARY <- "additive"
EK_PRIMARY        <- "Ek_min"
TOP_DECILE        <- 0.90
BOTTOM_DECILE     <- 0.10
METRICS           <- c("cyclomatic_complexity", "indeg", "btw")
EXPLORATORY_THRESHOLD <- 15L

joined_files <- list.files(out_dir, pattern = "_paths_joined\\.csv$", full.names = TRUE)
if (!length(joined_files)) stop("No joined files - run 01_join_static_dynamic.R first.")

paths_df <- load_full_paths_df()
node_df  <- load_full_node_df()

# Total statement count proxy: sum of cyclomatic complexity along the path
# (lacking a separate SLOC-per-path source). Documented in pre-reg deviation log
# if used differently.

matched_results <- list()
diag_smd <- list()

for (f in joined_files) {

  dt <- fread(f)
  m  <- unique(dt$model)
  dt <- dt[risk_form == RISK_FORM_PRIMARY & any_exercised == TRUE &
           !is.na(get(EK_PRIMARY))]

  if (!nrow(dt) || all(is.na(dt[[EK_PRIMARY]]))) {
    blocked_on_counts(paste0("matched contrast / model ", m))
    next
  }

  qhi <- quantile(dt$P_k_q50, TOP_DECILE,    na.rm = TRUE)
  qlo <- quantile(dt$P_k_q50, BOTTOM_DECILE, na.rm = TRUE)
  cv_med <- median(dt$CV_k, na.rm = TRUE)

  dt[, group := fifelse(P_k_q50 >= qhi & CV_k < cv_med, 1L,
                fifelse(P_k_q50 <= qlo, 0L, NA_integer_))]
  dt_match <- dt[!is.na(group)]

  if (sum(dt_match$group == 1L) < 5L || sum(dt_match$group == 0L) < 5L) {
    cat("[matched contrast] ", m, ": insufficient cases (",
        sum(dt_match$group == 1L), " high, ",
        sum(dt_match$group == 0L), " low). Skipping.\n", sep = "")
    next
  }

  # Path-mean static metrics for matched paths
  per_path_nodes <- expand_paths_to_nodes(m, RISK_FORM_PRIMARY,
                                          paths_df = paths_df,
                                          node_df  = node_df)
  metric_dt <- per_path_nodes[
    path_id %in% dt_match$path_id,
    lapply(.SD, mean, na.rm = TRUE),
    by = .(path_id),
    .SDcols = METRICS
  ]
  dt_match <- merge(dt_match, metric_dt, by = "path_id", all.x = TRUE)
  dt_match[, statement_proxy := cyclomatic_complexity * hops]

  # Mahalanobis NN matching with replacement (handover sec 3 Analysis 4)
  mout <- tryCatch(
    MatchIt::matchit(
      group ~ get(EK_PRIMARY) + hops + statement_proxy,
      data    = dt_match,
      method  = "nearest",
      distance = "mahalanobis",
      replace  = TRUE
    ),
    error = function(e) { warning("matchit failed for ", m, ": ", e$message); NULL }
  )
  if (is.null(mout)) next

  matched <- MatchIt::match.data(mout)
  n_pairs <- sum(matched$group == 1L)

  smd_pre  <- summary(mout)$sum.all[, "Std. Mean Diff."]
  smd_post <- summary(mout)$sum.matched[, "Std. Mean Diff."]
  diag_smd[[m]] <- data.table(
    model    = m,
    covariate = names(smd_post),
    smd_pre  = smd_pre[names(smd_post)],
    smd_post = smd_post
  )

  per_metric <- lapply(METRICS, function(mt) {
    a <- matched[matched$group == 1L, mt]
    b <- matched[matched$group == 0L, mt]
    if (!length(a) || !length(b)) return(NULL)
    # Use unpaired here since MatchIt with replace=TRUE breaks strict pairing.
    w <- suppressWarnings(wilcox.test(a, b))
    A <- effsize::VD.A(a, b)$estimate
    data.table(
      model    = m,
      metric   = mt,
      n_pairs  = n_pairs,
      median_high = median(a, na.rm = TRUE),
      median_low  = median(b, na.rm = TRUE),
      wilcox_p = w$p.value,
      effect_A = A,
      exploratory = n_pairs < EXPLORATORY_THRESHOLD
    )
  })
  matched_results[[m]] <- rbindlist(per_metric, fill = TRUE)
}

if (length(matched_results)) {
  fwrite(rbindlist(matched_results, fill = TRUE),
         file.path(syn_dir, "04_matched_contrast.csv"))
  fwrite(rbindlist(diag_smd, fill = TRUE),
         file.path(syn_dir, "04_matched_smd_diagnostics.csv"))
  print(rbindlist(matched_results, fill = TRUE))
  cat("\n[05_matched_contrast] wrote tables to ", syn_dir, "\n", sep = "")
} else {
  cat("\n[05_matched_contrast] no models had runtime counts - nothing written.\n")
}
