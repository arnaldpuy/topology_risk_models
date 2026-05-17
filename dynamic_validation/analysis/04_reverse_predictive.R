# 04_reverse_predictive.R - Analysis 3 (reverse predictive check)
# Holding execution frequency fixed, does static topology still discriminate?
# Among the top decile of paths by E_k, compare path-mean cyclomatic complexity,
# in-degree and betweenness between top- and bottom-P_k strata.
#
# BLOCKED ON RUNTIME COUNTS: requires Ek_min (or Ek_sum) for the runtime decile
# split. Also needs per-node metrics joined to the path - read from
# full_nodes_df.xlsx via path_id.

source("dynamic_validation/analysis/00_setup.R")

RISK_FORM_PRIMARY <- "additive"
EK_PRIMARY        <- "Ek_sum"  # see pre_registration.md deviation log 2026-05-17
TOP_RUNTIME_Q     <- 0.90
TOP_P_Q           <- 0.50  # within runtime-heavy paths, split at median P_k
METRICS           <- c("cyclomatic_complexity", "indeg", "btw")

joined_files <- list.files(out_dir, pattern = "_paths_joined\\.csv$", full.names = TRUE)
if (!length(joined_files)) stop("No joined files - run 01_join_static_dynamic.R first.")

# Per-path-per-node attributes built by joining full_paths_df.path_nodes to
# full_node_df. Cached per model below.
paths_df <- load_full_paths_df()
node_df  <- load_full_node_df()

reverse_results <- list()
plots <- list()

for (f in joined_files) {

  dt <- fread(f)
  m  <- unique(dt$model)
  dt <- dt[risk_form == RISK_FORM_PRIMARY & any_exercised == TRUE &
           !is.na(get(EK_PRIMARY))]

  if (!nrow(dt) || all(is.na(dt[[EK_PRIMARY]]))) {
    blocked_on_counts(paste0("reverse predictive / model ", m))
    next
  }

  # Top decile by runtime intensity
  q_runtime  <- quantile(dt[[EK_PRIMARY]], TOP_RUNTIME_Q, na.rm = TRUE)
  hot_paths  <- dt[get(EK_PRIMARY) >= q_runtime, .(model, path_id, P_k_q50)]

  # Path-mean static metrics built by expanding paths to (path_id, node) and
  # averaging the node-level attributes.
  per_path_nodes <- expand_paths_to_nodes(m, RISK_FORM_PRIMARY,
                                          paths_df = paths_df,
                                          node_df  = node_df)
  per_path_metrics <- per_path_nodes[
    path_id %in% hot_paths$path_id,
    lapply(.SD, mean, na.rm = TRUE),
    by = .(model, path_id),
    .SDcols = METRICS
  ]

  merged <- merge(hot_paths, per_path_metrics, by = c("model", "path_id"))
  if (nrow(merged) < 8) {
    cat("[reverse predictive] ", m, ": only ", nrow(merged),
        " runtime-heavy paths - skipping (need >=8).\n", sep = "")
    next
  }

  q_pk <- quantile(merged$P_k_q50, TOP_P_Q)
  merged[, stratum := ifelse(P_k_q50 >= q_pk, "high_Pk", "low_Pk")]

  per_metric <- lapply(METRICS, function(mt) {
    a <- merged[stratum == "high_Pk", get(mt)]
    b <- merged[stratum == "low_Pk",  get(mt)]
    if (!length(a) || !length(b)) return(NULL)
    w  <- suppressWarnings(wilcox.test(a, b))
    A  <- effsize::VD.A(a, b)$estimate
    data.table(model = m, metric = mt,
               n_high = length(a), n_low = length(b),
               median_high = median(a, na.rm = TRUE),
               median_low  = median(b, na.rm = TRUE),
               wilcox_p    = w$p.value,
               effect_A    = A)
  })
  reverse_results[[m]] <- rbindlist(per_metric, fill = TRUE)

  plot_dt <- melt(merged, id.vars = c("model", "path_id", "stratum"),
                  measure.vars = METRICS,
                  variable.name = "metric", value.name = "path_mean")
  plots[[m]] <- ggplot(plot_dt,
                       aes(x = stratum, y = path_mean, fill = stratum)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c("high_Pk" = "#C65D09", "low_Pk" = "#3B5B92")) +
    labs(title = m, x = NULL, y = "path-mean") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none")
}

if (length(reverse_results)) {
  out_tbl <- rbindlist(reverse_results, fill = TRUE)
  fwrite(out_tbl, file.path(syn_dir, "03_reverse_predictive.csv"))
  print(out_tbl)
  combined <- cowplot::plot_grid(plotlist = plots, ncol = 1, align = "v")
  ggsave(file.path(syn_dir, "03_reverse_predictive_boxplots.pdf"),
         combined, width = 8, height = 3 * length(plots))
  cat("\n[04_reverse_predictive] wrote table + boxplots to ", syn_dir, "\n", sep = "")
} else {
  cat("\n[04_reverse_predictive] no models had runtime counts - nothing written.\n")
}
