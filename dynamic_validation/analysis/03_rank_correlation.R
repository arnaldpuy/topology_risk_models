# 03_rank_correlation.R - Analysis 2 (Spearman P_k vs E_k)
# Among exercised paths, does the static risk score sort runtime intensity?
#
# BLOCKED ON RUNTIME COUNTS: needs Ek_min / Ek_sum to be non-NA, which requires
# Federico's runtime_callgraph_edges.csv to include a `runtime_calls` column.
# Until then the script reports skipped per model and writes no outputs.

source("dynamic_validation/analysis/00_setup.R")

RISK_FORM_PRIMARY <- "additive"
EK_PRIMARY        <- "Ek_min"   # min aggregator is the pre-registered primary
N_BOOT            <- 1000L
SEED              <- 123L

joined_files <- list.files(out_dir, pattern = "_paths_joined\\.csv$", full.names = TRUE)
if (!length(joined_files)) stop("No joined files - run 01_join_static_dynamic.R first.")

set.seed(SEED)

corr_results <- rbindlist(lapply(joined_files, function(f) {

  dt <- fread(f)
  m  <- unique(dt$model)
  dt <- dt[risk_form == RISK_FORM_PRIMARY & any_exercised == TRUE &
           !is.na(P_k_q50) & !is.na(get(EK_PRIMARY))]

  if (!nrow(dt) || all(is.na(dt[[EK_PRIMARY]]))) {
    blocked_on_counts(paste0("rank correlation / model ", m))
    return(NULL)
  }

  rho <- suppressWarnings(cor(dt$P_k_q50, dt[[EK_PRIMARY]], method = "spearman"))
  tau <- suppressWarnings(cor(dt$P_k_q50, dt[[EK_PRIMARY]], method = "kendall"))

  boot_rho <- boot::boot(
    data      = dt,
    statistic = function(d, i) suppressWarnings(
      cor(d$P_k_q50[i], d[[EK_PRIMARY]][i], method = "spearman")
    ),
    R = N_BOOT
  )
  ci <- tryCatch(
    boot::boot.ci(boot_rho, type = "perc")$percent[1, c(4, 5)],
    error = function(e) c(NA_real_, NA_real_)
  )

  # Diagnostic scatter plot
  pl <- ggplot(dt, aes(x = P_k_q50, y = log10(get(EK_PRIMARY) + 1))) +
    geom_point(alpha = 0.4, size = 0.8) +
    geom_smooth(method = "loess", se = FALSE, colour = "#C65D09") +
    labs(title  = sprintf("%s | Spearman rho = %.2f (95%% CI [%.2f, %.2f])",
                          m, rho, ci[1], ci[2]),
         x = expression(P[k]~"(median across ensemble)"),
         y = expression(log[10]*"("*E[k]+1*")")) +
    theme_minimal(base_size = 10)
  ggsave(file.path(syn_dir, paste0("02_scatter_", m, ".pdf")), pl,
         width = 5, height = 4)

  data.table(
    model     = m,
    n_paths   = nrow(dt),
    spearman  = rho,
    kendall   = tau,
    rho_ci_lo = ci[1],
    rho_ci_hi = ci[2]
  )
}), fill = TRUE)

if (nrow(corr_results)) {
  fwrite(corr_results, file.path(syn_dir, "02_rank_correlation.csv"))
  print(corr_results)
  cat("\n[03_rank_correlation] wrote table + per-model scatters to ", syn_dir, "\n", sep = "")
} else {
  cat("\n[03_rank_correlation] no models had runtime counts - nothing written.\n")
}
