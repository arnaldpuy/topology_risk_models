# Dynamic-validation pipeline

This folder implements the dynamic-profiling validation of the path-level
software-risk framework described in
[`../handover_dynamic_validation.md`](../handover_dynamic_validation.md). The
pre-committed decisions are locked in [`pre_registration.md`](pre_registration.md);
do not edit them without recording a deviation.

## Layout

```
dynamic_validation/
  ORCHIDEE/, PCR/, ...        # raw exports from the engineers (untracked)
  pre_registration.md         # locked decisions (sec 1-8 of handover)
  analysis/
    00_setup.R                # shared helpers, model registry, loaders
    01_join_static_dynamic.R  # build per-model paths_joined.csv
    02_coverage.R             # Analysis 1 - path coverage
    03_rank_correlation.R     # Analysis 2 - Spearman P_k vs E_k     [needs counts]
    04_reverse_predictive.R   # Analysis 3 - reverse predictive check [needs counts]
    05_matched_contrast.R     # Analysis 4 - matched-frequency       [needs counts]
    06_dormant_paths.R        # Analysis 5 - dormant high-risk paths
    07_cross_model_synthesis.R
  results/
    per_model/                # one *_paths_joined.csv per model
    synthesis/                # cross-model tables + main-text figures
```

## Inputs expected at repository root

These come from the existing static pipeline (`code_hydrological_models.R`):

- `full_paths_df.xlsx` - per-path point-estimate risk scores.
- `full_nodes_df.xlsx` - per-path nodes with their node-level metrics.
- `full_node_df.xlsx`  - per-model node attributes (file, function, C, indeg, btw).
- `full_ua_df.csv`     - per-path uncertainty summary (P_k_q50, P_k_q025/975 ...). Large (~1.7 GB); `00_setup.R` reads only the summary columns.
- `datasets/call_metrics/call_metrics_<lang>_<MODEL>.csv` - static edge lists.

And for each model with runtime data, in `dynamic_validation/<MODEL>/`:

- `runtime_callgraph_edges.csv` - runtime-observed edges. Currently lacks a `runtime_calls` column; once added, the per-path intensity `E_k` flows automatically and Analyses 2-4 become runnable.

## Running the pipeline

From the repository root:

```r
source("dynamic_validation/analysis/01_join_static_dynamic.R")
source("dynamic_validation/analysis/02_coverage.R")
source("dynamic_validation/analysis/06_dormant_paths.R")

# Once runtime_calls is shipped:
source("dynamic_validation/analysis/03_rank_correlation.R")
source("dynamic_validation/analysis/04_reverse_predictive.R")
source("dynamic_validation/analysis/05_matched_contrast.R")

source("dynamic_validation/analysis/07_cross_model_synthesis.R")
```

Each script writes to `dynamic_validation/results/`. Blocked analyses report
`[label] No per-edge runtime counts ... Skipping.` and exit cleanly.

## Adding a model

1. Drop the engineer's raw export under `dynamic_validation/<NEWMODEL>/` (must contain `runtime_callgraph_edges.csv` and `runtime_call_metrics_<lang>.csv`).
2. Append a row to `model_registry` in [`analysis/00_setup.R`](analysis/00_setup.R).
3. Re-run the pipeline.

## Current data status (as of pre-registration date)

| Model     | Runtime data shipped? | Per-edge counts? |
|-----------|-----------------------|------------------|
| ORCHIDEE  | yes (Fortran)         | no               |
| PCR-GLOBWB| yes (Python)          | no               |
| VIC       | not yet               | -                |
| HYPE      | not yet               | -                |

Engineer ask: please re-export `runtime_callgraph_edges.csv` with a
`runtime_calls` integer column giving the observed call count per
(caller, callee) edge. Schema otherwise unchanged.
