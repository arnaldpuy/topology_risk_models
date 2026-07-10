knitr::opts_chunk$set(echo = TRUE, dev = "pdf", cache = TRUE)


# PRELIMINARY FUNCTIONS ########################################################
################################################################################

sensobol::load_packages(c("data.table", "tidyverse", "openxlsx", "scales",
                          "cowplot", "readxl", "ggrepel", "tidytext", "here",
                          "tidygraph", "igraph", "foreach", "parallel", "ggraph",
                          "tools", "purrr", "sensobol", "benchmarkme", "softwareRisk",
                          "boot", "MatchIt", "effsize"))

# Set seed ---------------------------------------------------------------------

seed <- 123
set.seed(seed)


# PATHS ########################################################################

repo_root <- getwd()
cg_dir    <- file.path(repo_root, "dynamic_validation", "callgraph_csv", "Callgraphs")
out_dir   <- file.path(repo_root, "dynamic_validation", "results", "per_model")
syn_dir   <- file.path(repo_root, "dynamic_validation", "results", "synthesis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(syn_dir, recursive = TRUE, showWarnings = FALSE)

# MODEL REGISTRY ###############################################################
# static  = drift-corrected static call graph (callgraph_edges.csv)
# runtime = matched runtime trace with per-edge runtime_calls
# cc      = model key in datasets/cyclomatic_complexity_functions.csv
# role    = "primary" | "robustness"

model_registry <- data.table(
  model   = c("ORCHIDEE", "ORCHIDEE+IOIPSL", "PCR-GLOBWB"),
  static  = c(file.path(cg_dir, "13_ORCHIDEE-MAN-2025_ORIGINAL", "callgraph_edges.csv"),
              file.path(cg_dir, "13_ORCHIDEE-MAN-2025-AND IOIPSL", "callgraph_edges.csv"),
              file.path(cg_dir, "14_PCR-GLOBWB_model", "callgraph_edges.csv")),
  runtime = c(file.path(cg_dir, "13_ORCHIDEE-MAN-2025-AND IOIPSL", "runtime_callgraph_edges.csv"),
              file.path(cg_dir, "13_ORCHIDEE-MAN-2025-AND IOIPSL", "runtime_callgraph_edges.csv"),
              file.path(cg_dir, "14_PCR-GLOBWB_model", "runtime_callgraph_edges.csv")),
  cc      = c("ORCHIDEE", "ORCHIDEE", "PCR-GLOBWB"),
  role    = c("primary", "robustness", "primary")
)

available_models <- function() {
  reg <- copy(model_registry)
  reg[, available := file.exists(static) & file.exists(runtime)]
  reg[available == TRUE]
}

# NAME NORMALISATION ###########################################################
# Lowercase and strip a leading class qualifier (e.g. BmiPCRGlobWB.set_value ->
# set_value) so static and runtime symbols and the cyclomatic-complexity table
# join on a common key.

norm_fun_name <- function(x) {
  out <- tolower(trimws(as.character(x)))
  sub(".*\\.", "", out)
}


TOP_DECILE        <- 0.90
BOTTOM_DECILE     <- 0.10
EK_PRIMARY        <- "Ek_sum"   # see deviation log in pre_registration.md
EK_ROBUSTNESS     <- "Ek_min"
N_BOOT            <- 1000L
A_THRESHOLD       <- 0.55
MIN_PATHS_PER_BIN <- 20L
METRICS           <- c("pm_cc", "pm_indeg", "pm_btw")  # path-mean C / indeg / btw

# alpha/beta/gamma are the manuscript's node-risk weights.
ALPHA <- 0.6; BETA <- 0.3; GAMMA <- 0.1

color_strata <- c("top_decile" = "#C65D09", "bottom_decile" = "#3B5B92")


# BUILD STATIC P_k #############################################################
################################################################################

build_static <- function(edge_file, cc_model) {
  d  <- fread(edge_file)
  el <- unique(data.table(from = norm_fun_name(d$caller),
                          to   = norm_fun_name(d$callee)))
  el <- el[from != "" & to != "" & !is.na(from) & !is.na(to)]

  cc <- fread(file.path(repo_root, "datasets",
                        "cyclomatic_complexity_functions.csv"))[model == cc_model]
  cc[, name_norm := norm_fun_name(name)]
  cc_u    <- cc[, .(cyclomatic_complexity = max(cyclomatic_complexity, na.rm = TRUE)),
                by = name_norm]
  mean_cc <- mean(cc_u$cyclomatic_complexity, na.rm = TRUE)

  g   <- as_tbl_graph(el, directed = TRUE)
  nm  <- igraph::V(as.igraph(g))$name
  ccm <- cc_u[match(nm, name_norm), cyclomatic_complexity]
  ccm[is.na(ccm)] <- mean_cc
  g <- g |> activate(nodes) |> mutate(cyclomatic_complexity = ccm)

  out <- all_paths_fun(graph = g, alpha = ALPHA, beta = BETA, gamma = GAMMA,
                       p = 1, complexity_col = "cyclomatic_complexity")
  ua  <- uncertainty_fun(all_paths_out = out, N = 2^10, order = "first")

  paths <- as.data.table(out$paths)
  uap   <- as.data.table(ua$paths)
  uap[, P_k_q50  := vapply(uncertainty_analysis, median, numeric(1), na.rm = TRUE)]
  uap[, P_k_mean := vapply(uncertainty_analysis, mean,   numeric(1), na.rm = TRUE)]
  uap[, P_k_sd   := vapply(uncertainty_analysis, sd,     numeric(1), na.rm = TRUE)]
  uap[, CV_k     := P_k_sd / pmax(P_k_mean, 1e-12)]
  paths <- merge(paths, uap[, .(path_id, P_k_q50, P_k_mean, P_k_sd, CV_k)],
                 by = "path_id")

  list(paths = paths, nodes = as.data.table(out$nodes))
}


# JOIN RUNTIME E_k #############################################################
################################################################################

load_runtime <- function(runtime_file) {
  rt <- fread(runtime_file)
  rt[, `:=`(cn = norm_fun_name(caller), ce = norm_fun_name(callee))]
  if (!"runtime_calls" %in% names(rt)) rt[, runtime_calls := NA_integer_]
  rt <- rt[, .(runtime_calls = sum(runtime_calls, na.rm = TRUE)), by = .(cn, ce)]
  setkey(rt, cn, ce)
  rt
}

join_runtime <- function(static, runtime_file) {
  paths <- copy(static$paths); nodes <- static$nodes
  rt <- load_runtime(runtime_file)

  res <- paths[, {
    ns <- unlist(path_nodes)
    if (length(ns) < 2L) {
      list(n_edges = 0L, n_ex = 0L, frac_exercised = NA_real_,
           any_exercised = NA, all_exercised = NA,
           Ek_min = NA_real_, Ek_sum = NA_real_)
    } else {
      ed <- data.table(cn = norm_fun_name(ns[-length(ns)]),
                       ce = norm_fun_name(ns[-1L]))
      j  <- rt[ed, on = c("cn", "ce")]
      rc <- j$runtime_calls; ex <- !is.na(rc); rc[is.na(rc)] <- 0
      list(n_edges = nrow(ed), n_ex = sum(ex), frac_exercised = sum(ex) / nrow(ed),
           any_exercised = any(ex), all_exercised = all(ex),
           Ek_min = min(rc), Ek_sum = sum(rc))
    }
  }, by = .(path_id)]
  paths <- merge(paths, res, by = "path_id")

  pm <- paths[, {
    ns  <- norm_fun_name(unlist(path_nodes))
    sub <- nodes[name %in% ns]
    list(pm_cc    = mean(sub$cyclomatic_complexity, na.rm = TRUE),
         pm_indeg = mean(sub$indeg, na.rm = TRUE),
         pm_btw   = mean(sub$btw,   na.rm = TRUE))
  }, by = .(path_id)]
  merge(paths, pm, by = "path_id")
}

# Build + join for every available model.
reg <- available_models()
joined <- lapply(seq_len(nrow(reg)), function(i) {
  st <- build_static(reg$static[i], reg$cc[i])
  jn <- join_runtime(st, reg$runtime[i])
  jn[, model := reg$model[i]]
  fwrite(jn[, .(model, path_id, hops, P_k_q50, CV_k, n_edges, n_ex,
               frac_exercised, any_exercised, all_exercised, Ek_min, Ek_sum,
               pm_cc, pm_indeg, pm_btw)],
         file.path(out_dir, paste0(gsub("[^A-Za-z0-9]", "_", reg$model[i]),
                                   "_paths_joined.csv")))
  jn
})
names(joined) <- reg$model

edge_diag <- rbindlist(lapply(seq_len(nrow(reg)), function(i) {
  jn <- joined[[i]]
  data.table(model = reg$model[i], role = reg$role[i],
             paths = nrow(jn), exercised = sum(jn$any_exercised, na.rm = TRUE),
             median_hops = median(jn$hops, na.rm = TRUE))
}))
print(edge_diag)


# ANALYSIS 1 - PATH COVERAGE ###################################################
################################################################################

coverage_one <- function(dt, model) {
  dt <- dt[!is.na(P_k_q50)]
  qhi <- quantile(dt$P_k_q50, TOP_DECILE,    na.rm = TRUE)
  qlo <- quantile(dt$P_k_q50, BOTTOM_DECILE, na.rm = TRUE)
  cvm <- median(dt$CV_k, na.rm = TRUE)
  hi  <- dt$P_k_q50 >= qhi & dt$CV_k < cvm
  lo  <- dt$P_k_q50 <= qlo
  rbind(
    data.table(model = model, stratum = "top_decile", n = sum(hi),
               frac_any  = mean(dt$any_exercised[hi],  na.rm = TRUE),
               frac_edge = mean(dt$frac_exercised[hi], na.rm = TRUE),
               frac_all  = mean(dt$all_exercised[hi],  na.rm = TRUE),
               median_hops = as.numeric(median(dt$hops[hi], na.rm = TRUE))),
    data.table(model = model, stratum = "bottom_decile", n = sum(lo),
               frac_any  = mean(dt$any_exercised[lo],  na.rm = TRUE),
               frac_edge = mean(dt$frac_exercised[lo], na.rm = TRUE),
               frac_all  = mean(dt$all_exercised[lo],  na.rm = TRUE),
               median_hops = as.numeric(median(dt$hops[lo], na.rm = TRUE)))
  )
}

coverage_dt <- rbindlist(lapply(reg$model, function(m) coverage_one(joined[[m]], m)))
fwrite(coverage_dt, file.path(syn_dir, "01_coverage_table.csv"))
print(coverage_dt)


coverage_long <- melt(coverage_dt,
  id.vars = c("model", "stratum", "n", "median_hops"),
  measure.vars = c("frac_any", "frac_edge", "frac_all"),
  variable.name = "metric", value.name = "value")
coverage_long[, metric := factor(metric,
  levels = c("frac_any", "frac_edge", "frac_all"),
  labels = c("any edge", "mean fraction of edges", "all edges"))]

ggplot(coverage_long, aes(x = model, y = value, fill = stratum)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  facet_wrap(~ metric) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_manual(values = color_strata,
                    labels = c("bottom decile", "top decile")) +
  labs(x = NULL, y = "coverage", fill = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))


# ANALYSIS 1b - WITHIN-LENGTH DECILES ##########################################
################################################################################

within_length_one <- function(dt, model) {
  # Plain per-bin split-apply (avoids data.table by-group row recycling).
  x <- data.table(hops = dt$hops, P_k_q50 = dt$P_k_q50,
                  frac_exercised = dt$frac_exercised)[!is.na(P_k_q50)]
  x[, hb := cut(hops, c(0, 1, 2, 4, 7, Inf), labels = c("1", "2", "3-4", "5-7", "8+"))]
  rbindlist(lapply(levels(x$hb), function(b) {
    sub <- x[hb == b]
    if (nrow(sub) < MIN_PATHS_PER_BIN || length(unique(sub$P_k_q50)) < 2) return(NULL)
    q1 <- quantile(sub$P_k_q50, TOP_DECILE); q0 <- quantile(sub$P_k_q50, BOTTOM_DECILE)
    xt <- sub$frac_exercised[sub$P_k_q50 >= q1]; xb <- sub$frac_exercised[sub$P_k_q50 <= q0]
    if (length(xt) < 3 || length(xb) < 3) return(NULL)
    data.table(model = model, hb = b, n_top = length(xt), n_bot = length(xb),
      med_top = median(xt, na.rm = TRUE), med_bot = median(xb, na.rm = TRUE),
      wilcox_p = suppressWarnings(wilcox.test(xt, xb)$p.value),
      A = tryCatch(effsize::VD.A(xt, xb)$estimate, error = function(e) NA_real_))
  }))
}

within_dt <- rbindlist(lapply(reg$model, function(m) within_length_one(joined[[m]], m)),
                       fill = TRUE)
fwrite(within_dt, file.path(syn_dir, "01b_within_length.csv"))
print(within_dt)


wl_plot <- within_dt[!is.na(A)]
if (nrow(wl_plot)) {
  ggplot(wl_plot, aes(x = hb, y = A, fill = model)) +
    geom_col(position = position_dodge(0.7), width = 0.6) +
    geom_hline(yintercept = A_THRESHOLD, linetype = "dashed") +
    geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey60") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "path length (hops)", y = "effect size A (top vs bottom decile)",
         fill = NULL,
         caption = "dashed = A 0.55 threshold; dotted = A 0.5 (no discrimination)") +
    theme_minimal(base_size = 10) + theme(legend.position = "bottom")
}


# ANALYSIS 2 - RANK CORRELATION ################################################
################################################################################

rank_corr_one <- function(dt, model) {
  ex <- dt[any_exercised == TRUE & !is.na(get(EK_PRIMARY)) & !is.na(P_k_q50)]
  if (nrow(ex) < 10) return(NULL)
  rho  <- suppressWarnings(cor(ex$P_k_q50, ex[[EK_PRIMARY]],    method = "spearman"))
  tau  <- suppressWarnings(cor(ex$P_k_q50, ex[[EK_PRIMARY]],    method = "kendall"))
  rmin <- suppressWarnings(cor(ex$P_k_q50, ex[[EK_ROBUSTNESS]], method = "spearman"))
  bo <- boot::boot(ex, function(d, i)
    suppressWarnings(cor(d$P_k_q50[i], d[[EK_PRIMARY]][i], method = "spearman")),
    R = N_BOOT)
  ci <- tryCatch(boot::boot.ci(bo, type = "perc")$percent[1, c(4, 5)],
                 error = function(e) c(NA_real_, NA_real_))
  data.table(model = model, n_exercised = nrow(ex),
             spearman_sum = rho, ci_lo = ci[1], ci_hi = ci[2],
             kendall_sum = tau, spearman_min = rmin)
}

corr_dt <- rbindlist(lapply(reg$model, function(m) rank_corr_one(joined[[m]], m)),
                     fill = TRUE)
fwrite(corr_dt, file.path(syn_dir, "02_rank_correlation.csv"))
print(corr_dt)


scatter_panels <- lapply(reg$model, function(m) {
  ex <- joined[[m]][any_exercised == TRUE & !is.na(get(EK_PRIMARY))]
  if (!nrow(ex)) return(NULL)
  r <- corr_dt[model == m]
  ggplot(ex, aes(x = P_k_q50, y = log10(get(EK_PRIMARY) + 1))) +
    geom_point(alpha = 0.3, size = 0.7) +
    geom_smooth(method = "loess", se = FALSE, colour = "#C65D09") +
    labs(title = sprintf("%s (rho=%.2f [%.2f,%.2f], n=%d)",
                         m, r$spearman_sum, r$ci_lo, r$ci_hi, r$n_exercised),
         x = expression(P[k]), y = expression(log[10]*"("*E[k]^{sum}+1*")")) +
    theme_minimal(base_size = 9)
})
scatter_panels <- scatter_panels[!vapply(scatter_panels, is.null, logical(1))]
if (length(scatter_panels)) cowplot::plot_grid(plotlist = scatter_panels, nrow = 1)


# ANALYSIS 3 - REVERSE PREDICTIVE CHECK ########################################
################################################################################

reverse_one <- function(dt, model) {
  ex <- dt[any_exercised == TRUE & !is.na(get(EK_PRIMARY))]
  if (nrow(ex) < 16) return(NULL)
  qr  <- quantile(ex[[EK_PRIMARY]], TOP_DECILE, na.rm = TRUE)
  hot <- ex[get(EK_PRIMARY) >= qr]
  if (nrow(hot) < 8) return(NULL)
  qp <- quantile(hot$P_k_q50, 0.5)
  hot[, st := ifelse(P_k_q50 >= qp, "hi", "lo")]
  rbindlist(lapply(METRICS, function(mt) {
    a <- hot[st == "hi", get(mt)]; b <- hot[st == "lo", get(mt)]
    data.table(model = model, metric = mt, n_hi = length(a), n_lo = length(b),
               median_hi = median(a, na.rm = TRUE), median_lo = median(b, na.rm = TRUE),
               wilcox_p = suppressWarnings(wilcox.test(a, b)$p.value),
               A = tryCatch(effsize::VD.A(a, b)$estimate, error = function(e) NA_real_))
  }))
}

reverse_dt <- rbindlist(lapply(reg$model, function(m) reverse_one(joined[[m]], m)),
                        fill = TRUE)
reverse_dt[, metric := factor(metric, levels = METRICS,
                              labels = c("cyclomatic", "in-degree", "betweenness"))]
fwrite(reverse_dt, file.path(syn_dir, "03_reverse_predictive.csv"))
print(reverse_dt)


# ANALYSIS 4 - MATCHED-FREQUENCY CONTRAST ######################################
################################################################################

matched_one <- function(dt, model) {
  dt <- dt[!is.na(P_k_q50)]
  qhi <- quantile(dt$P_k_q50, TOP_DECILE,    na.rm = TRUE)
  qlo <- quantile(dt$P_k_q50, BOTTOM_DECILE, na.rm = TRUE)
  cvm <- median(dt$CV_k, na.rm = TRUE)
  dm  <- dt[any_exercised == TRUE & !is.na(get(EK_PRIMARY))]
  dm[, grp := fifelse(P_k_q50 >= qhi & CV_k < cvm, 1L,
              fifelse(P_k_q50 <= qlo,              0L, NA_integer_))]
  dm <- dm[!is.na(grp)]
  if (sum(dm$grp == 1) < 5 || sum(dm$grp == 0) < 5) return(NULL)
  dm[, stmt := pm_cc * hops]
  mo <- tryCatch(MatchIt::matchit(grp ~ Ek_sum + hops + stmt, data = dm,
                 method = "nearest", distance = "mahalanobis", replace = TRUE),
                 error = function(e) NULL)
  if (is.null(mo)) return(NULL)
  md <- MatchIt::match.data(mo); np <- sum(md$grp == 1)
  rbindlist(lapply(METRICS, function(mt) {
    a <- md[[mt]][md$grp == 1]; b <- md[[mt]][md$grp == 0]
    data.table(model = model, metric = mt, n_pairs = np,
               median_hi = median(a, na.rm = TRUE), median_lo = median(b, na.rm = TRUE),
               wilcox_p = suppressWarnings(wilcox.test(a, b)$p.value),
               A = tryCatch(effsize::VD.A(a, b)$estimate, error = function(e) NA_real_),
               exploratory = np < 15L)
  }))
}

matched_dt <- rbindlist(lapply(reg$model, function(m) matched_one(joined[[m]], m)),
                        fill = TRUE)
matched_dt[, metric := factor(metric, levels = METRICS,
                              labels = c("cyclomatic", "in-degree", "betweenness"))]
fwrite(matched_dt, file.path(syn_dir, "04_matched_contrast.csv"))
print(matched_dt)


# ANALYSIS 5 - DORMANT HIGH-RISK PATHS #########################################
################################################################################

dormant_one <- function(dt, model) {
  x <- data.table(hops = dt$hops, P_k_q50 = dt$P_k_q50, CV_k = dt$CV_k,
                  any_exercised = dt$any_exercised,
                  all_exercised = dt$all_exercised)[!is.na(P_k_q50)]
  x[, hb := cut(hops, c(0, 1, 2, 4, 7, Inf), labels = c("1", "2", "3-4", "5-7", "8+"))]
  rbindlist(lapply(levels(x$hb), function(b) {
    sub <- x[hb == b]
    if (nrow(sub) < MIN_PATHS_PER_BIN || length(unique(sub$P_k_q50)) < 2) return(NULL)
    q1 <- quantile(sub$P_k_q50, TOP_DECILE); cm <- median(sub$CV_k, na.rm = TRUE)
    h  <- sub[P_k_q50 >= q1 & CV_k < cm]
    data.table(model = model, hb = b, n_top = nrow(h),
               dormant_loose  = sum(h$any_exercised == FALSE | is.na(h$any_exercised)),
               dormant_strict = sum(h$all_exercised == FALSE | is.na(h$all_exercised)))
  }))
}

dormant_dt <- rbindlist(lapply(reg$model, function(m) dormant_one(joined[[m]], m)),
                        fill = TRUE)
dormant_dt[, `:=`(frac_loose  = dormant_loose  / n_top,
                  frac_strict = dormant_strict / n_top)]
fwrite(dormant_dt, file.path(syn_dir, "05_dormancy.csv"))
print(dormant_dt)


# CROSS-MODEL SYNTHESIS ########################################################
################################################################################

# Headline effect sizes on cyclomatic complexity from the two harder tests.
rev_cc <- reverse_dt[metric == "cyclomatic", .(model, A_reverse_cc = A)]
mat_cc <- matched_dt[metric == "cyclomatic", .(model, A_matched_cc = A, n_pairs)]
mat_btw <- matched_dt[metric == "betweenness", .(model, A_matched_btw = A)]

headline <- Reduce(function(a, b) merge(a, b, by = "model", all = TRUE),
                   list(corr_dt[, .(model, spearman_sum, ci_lo, ci_hi, n_exercised)],
                        rev_cc, mat_cc, mat_btw))
fwrite(headline, file.path(syn_dir, "06_cross_model_headline.csv"))
print(headline)

# Pre-registered decision rule: Spearman rho > 0 with CI_lo > 0, AND A > 0.55
# in at least one of the harder tests (reverse predictive or matched contrast).
decision <- headline[, .(model,
  spearman_pos = spearman_sum > 0 & ci_lo > 0,
  harder_A_pass = (A_reverse_cc > A_THRESHOLD) | (A_matched_cc > A_THRESHOLD))]
decision[, passes := spearman_pos & harder_A_pass]
fwrite(decision, file.path(syn_dir, "06_decision_rule.csv"))
print(decision)

n_primary <- reg[role == "primary", model]
cat(sprintf("\n[decision rule] %d of %d primary models pass both legs (%s).\n",
            sum(decision[model %in% n_primary, passes], na.rm = TRUE),
            length(n_primary), paste(n_primary, collapse = ", ")))
cat("Full pre-registered claim requires >= 3 of 4 target models (VIC + HYPE pending).\n")


p_rho <- ggplot(corr_dt, aes(x = model, y = spearman_sum)) +
  geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = expression("Spearman " ~ rho ~ "(P"[k]*", E"[k]^{sum}*")"),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 9) + theme(axis.text.x = element_text(angle = 20, hjust = 1))

effect_dt <- rbind(
  reverse_dt[metric == "cyclomatic", .(model, test = "reverse predictive", A)],
  matched_dt[metric == "cyclomatic", .(model, test = "matched contrast", A)])
p_A <- ggplot(effect_dt, aes(x = model, y = A, fill = test)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_hline(yintercept = A_THRESHOLD, linetype = "dashed") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Effect size A on cyclomatic complexity", x = NULL, y = NULL, fill = NULL) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))

cowplot::plot_grid(p_rho, p_A, ncol = 2, rel_widths = c(0.45, 0.55))


# SESSION INFORMATION ##########################################################
################################################################################

sessionInfo()

cat("Machine:     "); print(get_cpu()$model_name)
cat("Num cores:   "); print(detectCores(logical = FALSE))
cat("Num threads: "); print(detectCores(logical = TRUE))
