# 00_setup.R
# Shared helpers, model registry, and data loaders for the dynamic-validation
# pipeline. Sourced by every analysis script in this folder.
#
# Run from the repository root (the same convention as the rest of the repo).

sensobol::load_packages(c(
  "data.table", "tidyverse", "openxlsx", "readxl",
  "ggplot2", "cowplot", "scales",
  "boot", "MatchIt", "effsize"
))

# ---- Paths -------------------------------------------------------------------

repo_root <- getwd()
dv_dir    <- file.path(repo_root, "dynamic_validation")
raw_dir   <- function(model) file.path(dv_dir, model)
out_dir   <- file.path(dv_dir, "results", "per_model")
syn_dir   <- file.path(dv_dir, "results", "synthesis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(syn_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Model registry ----------------------------------------------------------
# `model` matches the keys used in `full_paths_df.xlsx` / `full_node_df.xlsx`.
# `raw_subdir` is the folder name under `dynamic_validation/`.
# `lang` records the language(s) of the model's runtime data; only models with
# at least one shipped runtime profile are listed.

model_registry <- data.table(
  model      = c("ORCHIDEE", "PCR-GLOBWB"),
  raw_subdir = c("ORCHIDEE", "PCR"),
  lang       = c("fortran",  "python")
)

# When VIC and HYPE arrive, append rows here. No code change required downstream.

available_models <- function() {
  reg <- copy(model_registry)
  reg[, available := vapply(raw_subdir, function(s) {
    f <- file.path(raw_dir(s), "runtime_callgraph_edges.csv")
    file.exists(f)
  }, logical(1))]
  reg[available == TRUE]
}

# ---- Static data loaders -----------------------------------------------------
# These read the canonical outputs produced by code_hydrological_models.R.
# Run that script (or restore the cached xlsx/csv) before running the dynamic
# validation pipeline.

load_full_paths_df <- function() {
  f <- file.path(repo_root, "full_paths_df.xlsx")
  if (!file.exists(f)) stop("Missing ", f,
    " - run code_hydrological_models.R first or sync from the main repo.")
  as.data.table(read.xlsx(f))
}

load_full_node_df <- function() {
  f <- file.path(repo_root, "full_node_df.xlsx")
  if (!file.exists(f)) stop("Missing ", f,
    " - run code_hydrological_models.R first or sync from the main repo.")
  as.data.table(read.xlsx(f))
}

load_full_nodes_df <- function() {
  f <- file.path(repo_root, "full_nodes_df.xlsx")
  if (!file.exists(f)) stop("Missing ", f,
    " - run code_hydrological_models.R first or sync from the main repo.")
  as.data.table(read.xlsx(f))
}

# Expand each static path into one row per (path_id, position, node) for a
# given model and risk_form, with node-level attributes (cyclomatic_complexity,
# indeg, btw, risk_score, file) attached. Used by analyses 3, 4, 5.

expand_paths_to_nodes <- function(model_label,
                                  risk_form_label = "additive",
                                  paths_df = NULL,
                                  node_df  = NULL) {
  if (is.null(paths_df)) paths_df <- load_full_paths_df()
  if (is.null(node_df))  node_df  <- load_full_node_df()

  fp <- paths_df[model == model_label & risk_form == risk_form_label]
  fn <- node_df[model == model_label,
                .(model, name, file, cyclomatic_complexity, indeg, btw,
                  risk_score)]

  if (!nrow(fp)) return(data.table())

  per_path <- rbindlist(lapply(seq_len(nrow(fp)), function(i) {
    nodes <- split_path_nodes(fp$path_nodes[i])[[1]]
    if (!length(nodes)) return(NULL)
    data.table(
      model   = model_label,
      path_id = fp$path_id[i],
      pos     = seq_along(nodes),
      name    = nodes
    )
  }), use.names = TRUE, fill = TRUE)

  merge(per_path, fn, by = c("model", "name"), all.x = TRUE)
}

# Per-path uncertainty summary (P_k_q50, P_k_mean, P_k_sd, CV_k, ...).
# Prefers the slim full_ua_summary.csv (produced by 00b_extract_ua_summary.R,
# ~tens of MB, with exact P_k_sd from the ensemble vector). Falls back to the
# raw full_ua_df.csv (~1.7 GB) with an sd approximated from quantiles.

load_uncertainty_summary <- function(models = NULL) {
  slim <- file.path(repo_root, "full_ua_summary.csv")
  raw  <- file.path(repo_root, "full_ua_df.csv")

  if (file.exists(slim)) {
    ua <- fread(slim)
  } else if (file.exists(raw)) {
    message("[setup] Slim full_ua_summary.csv not found; reading the raw ",
            "1.7 GB full_ua_df.csv. Run 00b_extract_ua_summary.R once to ",
            "speed this up.")
    cols <- c("model", "path_id", "path_nodes", "path_str", "hops",
              "P_k_min", "P_k_mean", "P_k_max",
              "P_k_q025", "P_k_q50", "P_k_q975")
    ua <- fread(raw, select = cols)
    ua[, P_k_sd := (P_k_q975 - P_k_q025) / (2 * 1.96)]
    ua[, CV_k   := P_k_sd / pmax(P_k_mean, 1e-12)]
  } else {
    stop("Missing both ", slim, " and ", raw,
         " - run code_hydrological_models.R UA section, then ",
         "00b_extract_ua_summary.R.")
  }

  if (!is.null(models)) ua <- ua[model %in% models]
  ua[]
}

# ---- Static edges per model --------------------------------------------------
# Reads the static call-metric CSVs and returns one normalised edge table per
# model with columns: caller_norm, callee_norm, model, lang.

load_static_edges <- function(model_label, lang) {
  fn <- sprintf("call_metrics_%s_%s.csv", lang, model_label)
  f  <- file.path(repo_root, "datasets", "call_metrics", fn)
  if (!file.exists(f)) stop("Missing ", f)
  dt <- fread(f)
  setnames(dt, c("function", "call"), c("caller", "callee"), skip_absent = TRUE)
  dt[, `:=`(
    caller_norm = norm_fun_name(caller),
    callee_norm = norm_fun_name(callee),
    model       = model_label,
    lang        = lang
  )]
  unique(dt[, .(caller_norm, callee_norm, model, lang)])
}

# ---- Runtime edge loader -----------------------------------------------------
# Reads runtime_callgraph_edges.csv for a given model. If a `runtime_calls`
# column is present (after Federico re-exports with counts), it is preserved;
# otherwise the loader returns a binary edge set with runtime_calls = NA.

load_runtime_edges <- function(model_label, raw_subdir) {
  f <- file.path(raw_dir(raw_subdir), "runtime_callgraph_edges.csv")
  if (!file.exists(f)) stop("Missing ", f)
  dt <- fread(f)
  # Schema as currently shipped: model, lang, caller, callee, caller_file,
  # caller_line, callee_file_hint, callee_kind. No runtime_calls column yet.
  dt[, `:=`(
    caller_norm = norm_fun_name(caller),
    callee_norm = norm_fun_name(callee)
  )]
  if (!"runtime_calls" %in% names(dt)) dt[, runtime_calls := NA_integer_]
  unique(dt[, .(caller_norm, callee_norm, runtime_calls)])
}

# ---- Name normalisation ------------------------------------------------------
# Both static and runtime function names need lowercasing and (for Python) class
# qualification stripped to maximise the joinable set. Adjust here if model-
# specific quirks emerge (Fortran name-mangling, Python decorators, etc.).

norm_fun_name <- function(x) {
  out <- tolower(trimws(as.character(x)))
  # Strip a leading class qualifier (e.g. "BmiPCRGlobWB.set_value" -> "set_value")
  out <- sub(".*\\.", "", out)
  out
}

# ---- Path edge derivation ----------------------------------------------------
# `full_paths_df.xlsx` stores path_nodes as a list-col that becomes a delimited
# string after the xlsx round-trip. Split it back into a character vector and
# derive the ordered edge list.

split_path_nodes <- function(path_nodes_str) {
  # Try a few common separators - the softwareRisk default is " -> ".
  for (sep in c(" -> ", ",", "|")) {
    if (any(grepl(sep, path_nodes_str, fixed = TRUE))) {
      return(strsplit(path_nodes_str, sep, fixed = TRUE))
    }
  }
  strsplit(path_nodes_str, "\\s+")
}

path_to_edges <- function(node_seq) {
  if (length(node_seq) < 2L) return(data.table(caller_norm = character(),
                                               callee_norm = character()))
  data.table(
    caller_norm = norm_fun_name(node_seq[-length(node_seq)]),
    callee_norm = norm_fun_name(node_seq[-1L])
  )
}

# ---- E_k aggregators ---------------------------------------------------------
# Given a per-path edge table joined to runtime counts (one row per edge of the
# path, with `runtime_calls` possibly NA for unmatched edges), compute the path-
# level intensity E_k under both aggregators.

aggregate_Ek <- function(edge_dt) {
  # Treat unmatched edges as 0 calls (cf. pre_registration.md sec 1).
  rc <- edge_dt$runtime_calls
  rc[is.na(rc)] <- 0
  if (length(rc) == 0L) return(list(Ek_min = NA_real_, Ek_sum = NA_real_))
  list(Ek_min = min(rc), Ek_sum = sum(rc))
}

# ---- Convenience -------------------------------------------------------------

has_runtime_counts <- function(runtime_edges) {
  any(!is.na(runtime_edges$runtime_calls))
}

blocked_on_counts <- function(label = "this analysis") {
  message("[", label, "] No per-edge runtime counts on the matched-static ",
          "edge set. Skipping. Re-run after the engineers re-export ",
          "runtime_callgraph_edges.csv with a `runtime_calls` column.")
  invisible(NULL)
}

cat("[setup] dynamic-validation helpers loaded. Models available:\n")
print(available_models())
