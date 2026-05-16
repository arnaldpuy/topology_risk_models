# 02_coverage.R - Analysis 1 (path coverage)
# Headline question: what fraction of statically-identified high-risk paths is
# exercised at runtime, vs the bottom-decile baseline?
#
# Reports three coverage metrics per stratum to defuse the path-length confound
# noted in the first pass (top-decile paths in ORCHIDEE are 5-10 hops; bottom
# decile are 1-4):
#   - frac_any  : at least one edge of the path is in the runtime trace
#   - frac_all  : every edge of the path is in the runtime trace (stricter)
#   - frac_edge : mean fraction of edges per path that are in the runtime trace
#     (length-insensitive)
#
# Also reports a length-stratified breakdown (coverage within hops bins) and
# runs both risk_forms (additive p=1, compensatory p=-1).
#
# Runs today on whatever models have a *_paths_joined.csv. Does not require
# runtime_calls.

source("dynamic_validation/analysis/00_setup.R")

# ---- Decision parameters (locked in pre_registration.md) ---------------------

TOP_DECILE    <- 0.90
BOTTOM_DECILE <- 0.10
RISK_FORMS    <- c("additive", "compensatory")

joined_files <- list.files(out_dir, pattern = "_paths_joined\\.csv$", full.names = TRUE)
if (!length(joined_files)) stop("No joined files - run 01_join_static_dynamic.R first.")

all_joined <- rbindlist(lapply(joined_files, fread), use.names = TRUE, fill = TRUE)

# ---- Headline coverage per model x risk_form x stratum -----------------------

coverage_dt <- all_joined[
  risk_form %in% RISK_FORMS & !is.na(P_k_q50),
  {
    qhi    <- quantile(P_k_q50, TOP_DECILE,    na.rm = TRUE)
    qlo    <- quantile(P_k_q50, BOTTOM_DECILE, na.rm = TRUE)
    cv_med <- median(CV_k, na.rm = TRUE)

    high_idx <- P_k_q50 >= qhi & CV_k < cv_med
    low_idx  <- P_k_q50 <= qlo

    rbind(
      data.table(
        stratum   = "top_decile",
        n_paths   = sum(high_idx),
        frac_any  = mean(any_exercised[high_idx],  na.rm = TRUE),
        frac_all  = mean(all_exercised[high_idx],  na.rm = TRUE),
        frac_edge = mean(frac_exercised[high_idx], na.rm = TRUE),
        median_hops = as.numeric(median(hops[high_idx], na.rm = TRUE))
      ),
      data.table(
        stratum   = "bottom_decile",
        n_paths   = sum(low_idx),
        frac_any  = mean(any_exercised[low_idx],  na.rm = TRUE),
        frac_all  = mean(all_exercised[low_idx],  na.rm = TRUE),
        frac_edge = mean(frac_exercised[low_idx], na.rm = TRUE),
        median_hops = as.numeric(median(hops[low_idx], na.rm = TRUE))
      )
    )
  },
  by = .(model, risk_form)
]

fwrite(coverage_dt, file.path(syn_dir, "01_coverage_table.csv"))

# ---- Length-stratified coverage (controls for path length) -------------------
# Bin paths by hops (1, 2, 3-4, 5-7, 8+) and report frac_any / frac_edge per
# (model, risk_form, hops_bin, stratum). When the same length bin appears in
# both top and bottom decile, the comparison is length-controlled.
#
# In ORCHIDEE + PCR-GLOBWB these bins did NOT overlap (P_k and hops are nearly
# co-linear under the saturating-OR aggregator). The within-length deciles
# block below (1b) does the apples-to-apples comparison within each hops bin
# directly and is the headline length-controlled result.

length_dt <- all_joined[
  risk_form %in% RISK_FORMS & !is.na(P_k_q50),
  {
    qhi    <- quantile(P_k_q50, TOP_DECILE,    na.rm = TRUE)
    qlo    <- quantile(P_k_q50, BOTTOM_DECILE, na.rm = TRUE)
    cv_med <- median(CV_k, na.rm = TRUE)

    .SD[, .(
      path_id,
      hops,
      any_exercised,
      frac_exercised,
      stratum = fifelse(P_k_q50 >= qhi & CV_k < cv_med, "top_decile",
                fifelse(P_k_q50 <= qlo,                 "bottom_decile",
                                                        NA_character_))
    )]
  },
  by = .(model, risk_form)
][!is.na(stratum)]

length_dt[, hops_bin := cut(hops,
                            breaks = c(0, 1, 2, 4, 7, Inf),
                            labels = c("1", "2", "3-4", "5-7", "8+"),
                            right  = TRUE)]

length_breakdown <- length_dt[
  ,
  .(n_paths   = .N,
    frac_any  = mean(any_exercised,  na.rm = TRUE),
    frac_edge = mean(frac_exercised, na.rm = TRUE)),
  by = .(model, risk_form, hops_bin, stratum)
][order(model, risk_form, hops_bin, stratum)]

fwrite(length_breakdown, file.path(syn_dir, "01_coverage_by_length.csv"))

# ---- Analysis 1b: within-length deciles --------------------------------------
# For each hops bin separately, take top vs bottom decile of P_k *within* the
# bin and compare coverage. This is the length-controlled headline metric.
#
# A bin is included only if it has >= MIN_PATHS_PER_BIN paths AND >= 2 distinct
# P_k_q50 values so the decile cut is meaningful.

MIN_PATHS_PER_BIN <- 20L

within_length_dt <- all_joined[risk_form %in% RISK_FORMS & !is.na(P_k_q50)]
within_length_dt[, hops_bin := cut(hops,
                                    breaks = c(0, 1, 2, 4, 7, Inf),
                                    labels = c("1", "2", "3-4", "5-7", "8+"),
                                    right  = TRUE)]

within_decile_dt <- within_length_dt[
  ,
  {
    if (.N < MIN_PATHS_PER_BIN || length(unique(P_k_q50)) < 2L) {
      data.table()
    } else {
      qhi <- quantile(P_k_q50, TOP_DECILE,    na.rm = TRUE)
      qlo <- quantile(P_k_q50, BOTTOM_DECILE, na.rm = TRUE)
      top <- P_k_q50 >= qhi
      bot <- P_k_q50 <= qlo
      rbind(
        data.table(
          stratum   = "top_decile",
          n_paths   = sum(top),
          frac_any  = mean(any_exercised[top],  na.rm = TRUE),
          frac_edge = mean(frac_exercised[top], na.rm = TRUE),
          frac_all  = mean(all_exercised[top],  na.rm = TRUE)
        ),
        data.table(
          stratum   = "bottom_decile",
          n_paths   = sum(bot),
          frac_any  = mean(any_exercised[bot],  na.rm = TRUE),
          frac_edge = mean(frac_exercised[bot], na.rm = TRUE),
          frac_all  = mean(all_exercised[bot],  na.rm = TRUE)
        )
      )
    }
  },
  by = .(model, risk_form, hops_bin)
]

fwrite(within_decile_dt, file.path(syn_dir, "01b_within_length_deciles.csv"))

# ---- Analysis 1b' : per-bin Wilcoxon + effect size A -------------------------
# For each (model, risk_form, hops_bin) cell with both strata present, do a
# top-vs-bottom Wilcoxon rank-sum on frac_exercised and report common-language
# effect size A. Captures within-length statistical discrimination.

within_test_dt <- within_length_dt[
  ,
  {
    if (.N < MIN_PATHS_PER_BIN || length(unique(P_k_q50)) < 2L) {
      data.table()
    } else {
      qhi <- quantile(P_k_q50, TOP_DECILE,    na.rm = TRUE)
      qlo <- quantile(P_k_q50, BOTTOM_DECILE, na.rm = TRUE)
      x_top <- frac_exercised[P_k_q50 >= qhi]
      x_bot <- frac_exercised[P_k_q50 <= qlo]
      if (length(x_top) < 3L || length(x_bot) < 3L) {
        data.table()
      } else {
        w <- suppressWarnings(wilcox.test(x_top, x_bot))
        A <- tryCatch(effsize::VD.A(x_top, x_bot)$estimate,
                      error = function(e) NA_real_)
        data.table(
          n_top      = length(x_top),
          n_bot      = length(x_bot),
          median_top = median(x_top, na.rm = TRUE),
          median_bot = median(x_bot, na.rm = TRUE),
          wilcox_p   = w$p.value,
          effect_A   = A
        )
      }
    }
  },
  by = .(model, risk_form, hops_bin)
]
within_test_dt[, signif := fcase(
  is.na(wilcox_p), "",
  wilcox_p < 0.001, "***",
  wilcox_p < 0.01,  "**",
  wilcox_p < 0.05,  "*",
  default = "ns"
)]

fwrite(within_test_dt, file.path(syn_dir, "01b_within_length_tests.csv"))

# ---- Plots -------------------------------------------------------------------

coverage_long <- melt(
  coverage_dt,
  id.vars       = c("model", "risk_form", "stratum", "n_paths", "median_hops"),
  measure.vars  = c("frac_any", "frac_all", "frac_edge"),
  variable.name = "metric",
  value.name    = "value"
)
coverage_long[, metric := factor(metric,
                                 levels = c("frac_any", "frac_edge", "frac_all"),
                                 labels = c("any edge exercised",
                                            "mean fraction of edges exercised",
                                            "all edges exercised"))]

p1 <- ggplot(coverage_long,
             aes(x = model, y = value, fill = stratum)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  facet_grid(risk_form ~ metric) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  scale_fill_manual(values = c("top_decile"    = "#C65D09",
                               "bottom_decile" = "#3B5B92"),
                    labels = c("top decile", "bottom decile")) +
  labs(x = NULL, y = "coverage", fill = NULL,
       title = "Analysis 1 - path coverage by static risk stratum",
       subtitle = "any vs mean-fraction-of-edges vs all (length-sensitive -> length-insensitive -> length-sensitive)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(syn_dir, "01_coverage.pdf"), p1, width = 9, height = 6)

p2 <- ggplot(length_breakdown,
             aes(x = hops_bin, y = frac_edge, fill = stratum)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  facet_grid(risk_form ~ model) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  scale_fill_manual(values = c("top_decile"    = "#C65D09",
                               "bottom_decile" = "#3B5B92"),
                    labels = c("top decile", "bottom decile")) +
  labs(x = "path length (hops)", y = "mean fraction of edges exercised",
       fill = NULL,
       title = "Analysis 1 - length-controlled coverage",
       subtitle = "shared hops bins give a length-controlled top-vs-bottom contrast") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(syn_dir, "01_coverage_by_length.pdf"), p2, width = 9, height = 6)

if (nrow(within_decile_dt)) {
  p3 <- ggplot(within_decile_dt,
               aes(x = hops_bin, y = frac_edge, fill = stratum)) +
    geom_col(position = position_dodge(0.7), width = 0.6) +
    facet_grid(risk_form ~ model) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1)) +
    scale_fill_manual(values = c("top_decile"    = "#C65D09",
                                 "bottom_decile" = "#3B5B92"),
                      labels = c("top decile", "bottom decile")) +
    labs(x = "path length (hops)", y = "mean fraction of edges exercised",
         fill = NULL,
         title = "Analysis 1b - within-length deciles of P_k",
         subtitle = "decile cuts taken *within* each hops bin (length is controlled by construction)") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
  ggsave(file.path(syn_dir, "01b_within_length_deciles.pdf"),
         p3, width = 9, height = 6)
} else {
  cat("[02_coverage] within-length deciles: no bin had >= ", MIN_PATHS_PER_BIN,
      " paths with variable P_k; skipping 1b plot.\n", sep = "")
}

# ---- Console summary ---------------------------------------------------------

cat("\n[02_coverage] headline (by model x risk_form x stratum):\n")
print(coverage_dt)
cat("\n[02_coverage] length-stratified (shared hops bins for top vs bottom):\n")
print(length_breakdown)
cat("\n[02_coverage] WITHIN-LENGTH deciles (top vs bottom of P_k inside each hops bin):\n")
print(within_decile_dt)
cat("\n[02_coverage] WITHIN-LENGTH Wilcoxon + effect size A on frac_exercised:\n")
print(within_test_dt)
cat("\n[02_coverage] wrote tables + plots to ", syn_dir, "\n", sep = "")
