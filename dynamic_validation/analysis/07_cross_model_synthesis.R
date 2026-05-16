# 07_cross_model_synthesis.R
# One summary table + main-text figure pulling the headline numbers from
# Analyses 1-3 across all available models.
#
# Reads the synthesis CSVs produced by 02-06. Gracefully degrades when files
# from blocked analyses (03-05) are missing.

source("dynamic_validation/analysis/00_setup.R")

safe_read <- function(path) {
  if (!file.exists(path)) {
    cat("[synth] missing ", path, " - skipping\n", sep = "")
    return(NULL)
  }
  fread(path)
}

coverage <- safe_read(file.path(syn_dir, "01_coverage_table.csv"))
corr     <- safe_read(file.path(syn_dir, "02_rank_correlation.csv"))
reverse  <- safe_read(file.path(syn_dir, "03_reverse_predictive.csv"))
matched  <- safe_read(file.path(syn_dir, "04_matched_contrast.csv"))
dormant  <- safe_read(file.path(syn_dir, "05_dormant_summary.csv"))

# ---- One-row-per-model headline table ----------------------------------------

models <- unique(c(coverage$model, corr$model, reverse$model,
                   matched$model, dormant$model))

headline <- data.table(model = models)

if (!is.null(coverage)) {
  cov_wide <- dcast(coverage, model ~ stratum,
                    value.var = c("frac_any", "frac_all"))
  headline <- merge(headline, cov_wide, by = "model", all.x = TRUE)
}

if (!is.null(corr)) {
  headline <- merge(headline,
                    corr[, .(model,
                             spearman_Pk_Ek    = spearman,
                             spearman_ci_lo    = rho_ci_lo,
                             spearman_ci_hi    = rho_ci_hi)],
                    by = "model", all.x = TRUE)
}

if (!is.null(reverse)) {
  rev_wide <- dcast(reverse, model ~ metric, value.var = "effect_A")
  setnames(rev_wide,
           setdiff(names(rev_wide), "model"),
           paste0("A_reverse_", setdiff(names(rev_wide), "model")))
  headline <- merge(headline, rev_wide, by = "model", all.x = TRUE)
}

if (!is.null(matched)) {
  m_wide <- dcast(matched, model ~ metric, value.var = "effect_A")
  setnames(m_wide,
           setdiff(names(m_wide), "model"),
           paste0("A_matched_", setdiff(names(m_wide), "model")))
  headline <- merge(headline, m_wide, by = "model", all.x = TRUE)
}

if (!is.null(dormant)) {
  headline <- merge(headline,
                    dormant[, .(model, frac_dormant_loose, frac_dormant_strict)],
                    by = "model", all.x = TRUE)
}

fwrite(headline, file.path(syn_dir, "06_cross_model_headline.csv"))
print(headline)

# ---- Pre-registered decision rule (pre_registration.md sec 6) ----------------

dec <- NULL
if (!is.null(corr) && !is.null(reverse)) {
  cor_ok <- corr[, .(spearman_pos = spearman > 0 & rho_ci_lo > 0), by = model]
  rev_ok <- reverse[, .(any_A_above = any(effect_A > 0.55, na.rm = TRUE)),
                   by = model]
  dec <- merge(cor_ok, rev_ok, by = "model", all = TRUE)
  dec[, supports_framework := spearman_pos & any_A_above]
  fwrite(dec, file.path(syn_dir, "06_decision_rule.csv"))
  cat("\n[synth] decision-rule per model:\n")
  print(dec)
  n_support <- sum(dec$supports_framework, na.rm = TRUE)
  cat("[synth] ", n_support, " of ", nrow(dec),
      " models meet the joint criterion ",
      "(>=3 of 4 required for full support).\n", sep = "")
}

# ---- Headline figure (placeholder layout) ------------------------------------

if (!is.null(coverage) && nrow(coverage)) {
  p_cov <- ggplot(coverage,
                  aes(x = model, y = frac_any, fill = stratum)) +
    geom_col(position = position_dodge(0.7), width = 0.6) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c("top_decile" = "#C65D09",
                                 "bottom_decile" = "#3B5B92")) +
    labs(title = "Coverage", x = NULL, y = NULL, fill = NULL) +
    theme_minimal(base_size = 10) + theme(legend.position = "bottom")

  panels <- list(p_cov)
  if (!is.null(corr) && nrow(corr)) {
    p_corr <- ggplot(corr, aes(x = model, y = spearman)) +
      geom_pointrange(aes(ymin = rho_ci_lo, ymax = rho_ci_hi)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = expression("Spearman " ~ rho ~ "(P"[k] ~ ", E"[k] ~ ")"),
           x = NULL, y = NULL) +
      theme_minimal(base_size = 10)
    panels <- c(panels, list(p_corr))
  }
  if (!is.null(reverse) && nrow(reverse)) {
    p_rev <- ggplot(reverse, aes(x = model, y = effect_A, fill = metric)) +
      geom_col(position = position_dodge(0.7), width = 0.6) +
      geom_hline(yintercept = 0.5, linetype = "dashed") +
      labs(title = "Reverse predictive (effect size A)",
           x = NULL, y = NULL, fill = NULL) +
      theme_minimal(base_size = 10) + theme(legend.position = "bottom")
    panels <- c(panels, list(p_rev))
  }

  combined <- cowplot::plot_grid(plotlist = panels, ncol = 1, align = "v")
  ggsave(file.path(syn_dir, "06_main_figure.pdf"), combined,
         width = 6, height = 3 * length(panels))
  cat("[synth] main figure written.\n")
}

cat("\n[07_cross_model_synthesis] done.\n")
