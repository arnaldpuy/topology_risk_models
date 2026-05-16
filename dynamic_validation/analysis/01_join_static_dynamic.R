# 01_join_static_dynamic.R
# For each model with runtime data, build a per-path joined table containing:
#   - static path identifiers and risk summary (P_k_q50, CV_k, hops, ...),
#   - per-path counts of edges exercised at runtime under the looser ("any")
#     and stricter ("all") definitions from pre_registration.md sec 2,
#   - per-path E_k under both aggregators (min, sum). NA until runtime counts
#     for the matched-static edge set arrive.
#
# Output: dynamic_validation/results/per_model/<model>_paths_joined.csv
#
# Run from repository root.

source("dynamic_validation/analysis/00_setup.R")

# ---- Load static side once ---------------------------------------------------

ua_summary  <- load_uncertainty_summary(models = available_models()$model)
node_df     <- load_full_node_df()

# Identify the risk_form column if present (the canonical file rbinds additive +
# compensatory). If absent, we assume additive.
if (!"risk_form" %in% names(ua_summary)) ua_summary[, risk_form := "additive"]

# ---- Per-model join ----------------------------------------------------------

join_one_model <- function(model_label, raw_subdir, lang) {

  cat("\n[join] model = ", model_label, "  lang = ", lang, "\n", sep = "")

  static_edges  <- load_static_edges(model_label, lang)
  runtime_edges <- load_runtime_edges(model_label, raw_subdir)

  # Matched-static runtime edge set: edges in the static call graph that also
  # appear in the runtime trace. This is the set whose counts we want from the
  # profiler. Currently runtime_calls is NA on this set.
  matched_edges <- merge(static_edges, runtime_edges,
                         by = c("caller_norm", "callee_norm"),
                         all.x = TRUE)
  matched_edges[, exercised := !is.na(runtime_calls) |
                  paste(caller_norm, callee_norm) %in%
                    paste(runtime_edges$caller_norm, runtime_edges$callee_norm)]

  # For each static path, derive its edge list and join the per-edge runtime
  # status / count.
  paths <- ua_summary[model == model_label]
  if (!nrow(paths)) {
    warning("No static paths found for ", model_label,
            " in full_ua_df.csv (or risk_form filter dropped them).")
    return(invisible(NULL))
  }

  cat("[join] ", nrow(paths), " paths in static; ",
      nrow(static_edges), " static edges; ",
      nrow(runtime_edges), " runtime edges; ",
      sum(matched_edges$exercised), " matched.\n", sep = "")

  # Build a fast lookup keyed by (caller_norm, callee_norm).
  setkey(matched_edges, caller_norm, callee_norm)

  per_path <- paths[, {
    node_seq <- split_path_nodes(path_nodes)[[1]]
    edges    <- path_to_edges(node_seq)
    if (nrow(edges) == 0L) {
      list(
        n_edges          = 0L,
        n_edges_static   = 0L,
        n_edges_exercised = 0L,
        frac_exercised   = NA_real_,
        any_exercised    = NA,
        all_exercised    = NA,
        Ek_min           = NA_real_,
        Ek_sum           = NA_real_
      )
    } else {
      joined <- matched_edges[edges, on = c("caller_norm", "callee_norm")]
      ex     <- !is.na(joined$exercised) & joined$exercised
      Ek     <- aggregate_Ek(joined)
      list(
        n_edges           = nrow(edges),
        n_edges_static    = sum(!is.na(joined$model)),
        n_edges_exercised = sum(ex),
        frac_exercised    = sum(ex) / nrow(edges),
        any_exercised     = any(ex),
        all_exercised     = all(ex),
        Ek_min            = Ek$Ek_min,
        Ek_sum            = Ek$Ek_sum
      )
    }
  }, by = .(model, path_id, risk_form, hops,
            P_k_min, P_k_mean, P_k_max, P_k_q025, P_k_q50, P_k_q975, CV_k)]

  out_f <- file.path(out_dir, paste0(model_label, "_paths_joined.csv"))
  fwrite(per_path, out_f)
  cat("[join] wrote ", out_f, "\n", sep = "")

  # Also persist the matched-edge diagnostics so downstream scripts and the
  # README can quote them.
  diag_f <- file.path(out_dir, paste0(model_label, "_edge_diagnostics.csv"))
  fwrite(data.table(
    model            = model_label,
    lang             = lang,
    static_edges     = nrow(static_edges),
    runtime_edges    = nrow(runtime_edges),
    matched_edges    = sum(matched_edges$exercised),
    runtime_coverage_pct = 100 * sum(matched_edges$exercised) / nrow(static_edges),
    has_runtime_counts = has_runtime_counts(runtime_edges)
  ), diag_f)

  invisible(per_path)
}

# ---- Run over every available model ------------------------------------------

reg <- available_models()
for (i in seq_len(nrow(reg))) {
  join_one_model(reg$model[i], reg$raw_subdir[i], reg$lang[i])
}

cat("\n[join] done. Joined tables in ", out_dir, "\n", sep = "")
