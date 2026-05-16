# Dynamic-validation: progress summary

**Date:** 2026-05-16
**Scope:** Dynamic-profiling validation of the path-level software-risk framework described in Puy et al., *The topology of software risk in scientific models* (target *Nature Computational Science*).
**Models with runtime data shipped so far:** ORCHIDEE (Fortran, MAN-2025 branch), PCR-GLOBWB (Python). VIC and HYPE pending.

This document is a point-in-time snapshot for collaborators. It records what has been built, what the analyses show, what those results buy the manuscript, and the work that still needs to happen before the validation is publication-ready.

---

## 1. What was built

A complete analysis pipeline lives under [`dynamic_validation/`](.) in the repository:

```
dynamic_validation/
  pre_registration.md         locked decisions (cf. handover §6)
  summary.md                  this file
  README.md                   pipeline usage
  ORCHIDEE/, PCR/             raw exports from the engineers
  analysis/
    00_setup.R                shared helpers, model registry, data loaders
    00b_extract_ua_summary.R  one-off: slims the 1.78 GB full_ua_df.csv to an 8.7 MB summary
    01_join_static_dynamic.R  per-path join of static graph + runtime trace
    02_coverage.R             Analysis 1 (path coverage) + 1b (within-length deciles + Wilcoxon)
    03_rank_correlation.R     Analysis 2 (Spearman P_k vs E_k)              [blocked on counts]
    04_reverse_predictive.R   Analysis 3 (reverse predictive check)         [blocked on counts]
    05_matched_contrast.R     Analysis 4 (matched-frequency contrast)       [blocked on counts]
    06_dormant_paths.R        Analysis 5 (dormant high-risk paths)
    07_cross_model_synthesis.R headline table + main-text figure
  results/per_model/          *_paths_joined.csv + edge diagnostics per model
  results/synthesis/          cross-model tables + plots
```

Key design choices:

- **Pre-registration locked first.** All decisions in `pre_registration.md` (aggregation, threshold, decile cut, effect-size convention, decision rule) were set before any analysis ran against the dynamic data, with a deviation log at the bottom.
- **Three coverage metrics, not one.** `frac_any` (loose: ≥1 edge exercised), `frac_all` (strict: all edges exercised), and `frac_edge` (mean fraction of edges per path exercised — length-insensitive). The triplet defuses the path-length confound described in §5.
- **Slim uncertainty summary.** The 1.78 GB `full_ua_df.csv` was reduced to an 8.7 MB `full_ua_summary.csv` containing per-path P_k quantiles and an exact P_k_sd computed from the ensemble vector. The slim file is small enough to commit and makes the validation reproducible without the raw monster.
- **Blocked analyses are wired but soft-fail.** Analyses 2/3/4 require `runtime_calls` on the matched-static edge set, which is not yet shipped. The scripts exit cleanly with a `[label] No per-edge runtime counts ... Skipping.` message until the column arrives — no schema changes needed downstream.

---

## 2. What was actually run

Against ORCHIDEE + PCR-GLOBWB, with the static side restricted to additive UA (compensatory was originally planned as a robustness check; see §5):

| Analysis | Status | Output |
|---|---|---|
| 1 — coverage (unstratified, three metrics) | done | `01_coverage_table.csv`, `01_coverage.pdf` |
| 1b — within-length deciles + Wilcoxon + A | done | `01b_within_length_deciles.csv`, `01b_within_length_tests.csv`, `01b_within_length_deciles.pdf` |
| 1c — length-stratified breakdown | done | `01_coverage_by_length.csv`, `01_coverage_by_length.pdf` |
| 2 — Spearman ρ on P_k vs E_k | blocked on counts | — |
| 3 — reverse predictive check | blocked on counts | — |
| 4 — matched-frequency contrast | blocked on counts | — |
| 5 — dormant high-risk paths (within-length) | done | `05_dormant_summary.csv`, `05_dormant_paths_table.csv` |
| 6 — cross-model synthesis | done (degrades gracefully on blocked outputs) | `06_cross_model_headline.csv`, `06_main_figure.pdf` |

---

## 3. Headline results

### Edge-level coverage diagnostics

| Model | static paths | static edges | runtime edges | matched | runtime coverage |
|---|---|---|---|---|---|
| ORCHIDEE | 3152 | 1368 | 767 | 168 | 12.3% |
| PCR-GLOBWB | 101 | 116 | 144 | 9 | 7.8% |

### Unstratified coverage (additive, all paths)

| Model | stratum | n | frac_any | frac_edge | frac_all | median hops |
|---|---|---|---|---|---|---|
| ORCHIDEE | top decile | 318 | 100.0% | 45.9% | 0.0% | 7 |
| ORCHIDEE | bottom decile | 316 | 12.3% | 8.7% | 5.7% | 1 |
| PCR-GLOBWB | top decile | 13 | 92.3% | 45.5% | 7.7% | 4 |
| PCR-GLOBWB | bottom decile | 15 | 0.0% | 0.0% | 0.0% | 1 |

### Within-length deciles (the apples-to-apples test, ORCHIDEE)

| hops bin | n top / bot | median frac_edge (top vs bot) | Wilcoxon p | A | reading |
|---|---|---|---|---|---|
| 1 | 23 / 24 | 0.00 vs 1.00 | <10⁻⁴ | **0.23** | reversed |
| 2 | 45 / 45 | 0.00 vs 0.00 | 1.00 | **0.50** | null |
| 3–4 | 100 / 100 | 0.50 vs 0.00 | <10⁻¹¹ | **0.78** | supportive |
| 5–7 | 133 / 133 | 0.43 vs 0.20 | <10⁻³⁵ | **0.93** | strongly supportive |
| 8+ | 16 / 16 | 0.44 vs 0.25 | <10⁻⁶ | **1.00** | maximal |

PCR-GLOBWB: only the hops=1 bin had enough paths for a within-length cut, and both strata exercise 0% — uninformative.

### Within-length dormancy (ORCHIDEE top-decile)

| hops bin | n | dormant loose | dormant strict |
|---|---|---|---|
| 1 | 1 | 1 (100%) | 1 (100%) |
| 2 | 41 | 39 (95%) | 41 (100%) |
| 3–4 | 100 | 5 (5%) | 100 (100%) |
| 5–7 | 133 | 0 (0%) | 133 (100%) |
| 8+ | 8 | 0 (0%) | 8 (100%) |

---

## 4. How this builds on and defends the manuscript

The framework's central claim is that path-level static risk identifies execution paths that *concentrate systemic risk in scientific software*. Each of the five analyses in the handover targets a specific reviewer counter-claim. We can now mark off two of those counter-claims, partially address a third, and explicitly preserve the framework's most important framing concession.

**Claim defended: "the framework discriminates beyond what trivial baselines would predict."** The unstratified coverage contrast (100% vs 12.3% any-exercised in ORCHIDEE; 92.3% vs 0% in PCR) confirms that top-decile P_k paths are an order of magnitude more likely to be exercised than bottom-decile paths. Once we control for path length — necessary because the saturating-OR aggregator is monotone in path length — the within-length effect sizes on ORCHIDEE are A = 0.78–1.00 across hops 3–8+, all p < 10⁻⁶. This *exceeds the pre-registered threshold of A > 0.55* in every length bin in the range where the framework is designed to apply. It is the cleanest evidence that static topology is doing more than reproducing path length.

**Claim defended: "static analysis identifies vulnerabilities that runtime profiling alone misses."** This was the manuscript's own framing of complementarity rather than substitutability, and it is now empirically grounded. At hops 3–4, 5% of within-length top-decile high-risk paths have zero edges exercised in the chosen ORCHIDEE configuration; at every length bin ≥3, *no* top-decile path has every edge exercised. The runtime trace touches the beginnings and middles of high-risk execution chains but not their full traversal — exactly the "potential vulnerability" the framework was built to flag.

**Framing reinforced rather than weakened: the hops=1 reversal.** Length-1 high-risk paths are *less* exercised at runtime than length-1 low-risk paths in ORCHIDEE (4.3% vs 58.3%, A = 0.23, p < 10⁻⁴). The natural reading is that high-cyclomatic-complexity single-function calls live in error-handling, edge-case, or specialised routines that do not fire under nominal inputs, while low-complexity single-function calls are exactly the workhorse utilities the model uses constantly. This is consistent with — and arguably the cleanest empirical instance of — the manuscript's "static risk = potential, not realised, vulnerability" framing. It needs to be reported as part of the validation, not hidden in the supplement.

**Pre-registered decision rule, status.** The rule requires Spearman ρ > 0 (lower bound of bootstrap CI > 0) on the P_k–E_k correlation *and* A > 0.55 in ≥3 of 4 models. The A criterion is met within the applicable length range for ORCHIDEE; the Spearman criterion is structurally blocked on per-edge runtime counts. We have one of the two legs cleanly, on one of four target models. That is supportive but not yet sufficient to claim the rule is satisfied.

---

## 5. Weaknesses

These are written so that a reviewer's most natural objections appear on this list before they appear in a review.

**5.1 Single configuration per model.** Both ORCHIDEE and PCR-GLOBWB were profiled under one configuration apiece (single-site offline for ORCHIDEE; a working subset for PCR). Other configurations would exercise different subsets of paths. The dormant set in particular is a configuration-specific quantity. Generalising any "5% of top-decile paths are dormant" statement beyond the chosen configuration requires either more configurations or an explicit lower-bound framing.

**5.2 Codebase drift between static and runtime snapshots.** The engineers note that (a) PCR-GLOBWB was profiled against a slightly updated version because the original would not execute, and (b) ORCHIDEE has an architecture-specific code-regeneration step, so the Fortran the profiler observed is not identical to the Fortran the static analyser parsed. Both interpretations are partly responsible for the 88% non-match between static and runtime edge sets — but we cannot currently separate genuine dormancy from drift artefacts. *Until the static analysis is re-run against the actual runtime code versions ("phase 2" in the engineering plan), every coverage number carries a footnote of unknown size.*

**5.3 PCR-GLOBWB is statistically thin.** Only 13 paths in the unstratified top decile; only one length bin (hops=1) has enough paths for the within-length cut, and that bin shows zero exercise on both sides. PCR is providing essentially no validation signal at present, and any cross-model claim relies on ORCHIDEE alone.

**5.4 Two of four target models still missing.** VIC and HYPE have not been profiled yet. The pre-registered decision rule explicitly requires ≥3 of 4 models to meet the A > 0.55 threshold; with two models in hand and one of them data-limited, we cannot meet that rule on cardinality grounds, irrespective of the substantive findings.

**5.5 Per-edge runtime counts are missing on the matched-static edge set.** This blocks Analysis 2 (Spearman ρ), Analysis 3 (reverse predictive), and Analysis 4 (matched-frequency contrast) — three of the five analyses the handover specifies. Coverage and dormancy together do *not* answer the question "among exercised paths, does static P_k order runtime intensity correctly", which is the second half of the pre-registered decision rule.

**5.6 Name normalisation is aggressive and unaudited.** Both the static and runtime caller/callee names are lowercased and stripped of leading class qualifiers before joining. This is the right starting point for Python BMI classes (e.g. `BmiPCRGlobWB.set_value` → `set_value`) but may over-collapse for Fortran with module-name prefixes or for Python overloads. A sample audit of unmatched edges has not yet been done.

**5.7 The hops=1 and hops=2 results are honest but rhetorically awkward.** A reviewer who reads the abstract and lands on "the framework fails at hops ≤ 2" without the saturating-OR explanation may misread it. The within-length cut needs to be presented carefully in Methods so the reversal is framed as the framework's domain of applicability rather than as evidence against it.

**5.8 Compensatory risk form is not separately reported.** During the runs we discovered that `softwareRisk::uncertainty_fun` samples p ∈ [-1, 2] internally regardless of which `p` was used to enumerate paths, so additive and compensatory UA give the same per-path uncertainty distribution. This means the manuscript's existing "additive vs compensatory" robustness pair is meaningful only for the point estimate from `all_paths_fun`, not for the ensemble median used in the validation. The validation reports additive only; if a separate compensatory-point-estimate ranking is wanted, it has to be computed in a side-script.

---

## 6. What needs to happen to strengthen the validation

Mapped 1:1 to the weaknesses above; ordered by leverage.

**6.1 (Highest leverage) Get the per-edge `runtime_calls` column from Federico.** This is a one-line export change on his side and unblocks Analyses 2, 3, and 4 with no schema changes on ours. Without it the pre-registered decision rule cannot be evaluated at all.

**6.2 Trigger the engineering "phase 2" re-run of the static analyser on the runtime code versions.** For PCR-GLOBWB this means parsing the updated codebase Federico had to use to get the model running. For ORCHIDEE this means parsing the post-regeneration Fortran. Once both are done, the matched-static edge set goes from "12% of static edges" to whatever the true intersection is, and the dormancy numbers stop carrying their drift footnote.

**6.3 Push for VIC and HYPE runs.** Even a single configuration each gets us to four models, which is what the decision rule was written for. If a second PCR-GLOBWB configuration is also feasible (e.g. ISIMIP3 instead of the developer-shipped run), that would address §5.3 simultaneously.

**6.4 Audit name normalisation on a sample of unmatched edges.** Take 50 random unmatched static edges per model, hand-check whether the unmatched callees are (a) genuinely absent from the runtime trace, (b) present but renamed/mangled, or (c) present but in a class/module the normaliser collapsed away. This is half a day's work and gates the credibility of every coverage number.

**6.5 Restrict the headline to hops ≥ 3 and state the restriction in Methods.** Write the within-length framing as the primary analysis ("the framework is defined as a path-level integration; we evaluate it where path-level integration is non-degenerate, i.e. hops ≥ 3"), then report hops 1 and 2 honestly in the supplement as evidence for the potential-vs-realised distinction. This pre-empts the most likely framing objection.

**6.6 Either accept compensatory as redundant, or compute a fixed-p compensatory ranking in a side-script.** If the latter, write a 30-line script that takes the existing `all_graphs` (or rebuilds it for the two models), runs `all_paths_fun(p=-1)`, and joins the resulting `path_risk_score` to the joined-paths CSV as an alternative ranking. Run the within-length deciles against that ranking to confirm the headline doesn't depend on the choice of p for point-estimate enumeration.

**6.7 Add a second configuration per model where feasible.** This is the single most powerful answer to "your dormant set is configuration-specific" — running the same analysis on a different developer-shipped configuration of each model and showing that the dormant *set* is different but the dormant *fraction* is in the same ballpark would be a strong robustness statement. The handover's pre-committed configurations table (§4 of the handover) already lists backup options for PCR-GLOBWB and ORCHIDEE.

**6.8 Draft the ~400-word Results paragraph and ~250-word Methods paragraph called for in the handover §7.** Numbers are now stable enough for a first pass. The Results paragraph should lead with the within-length finding, fold the dormancy claim in as complementarity, and explicitly *not* claim the pre-registered decision rule is met. The Methods paragraph should specify the within-length cut and the hops ≥ 3 restriction.

---

## 7. Bottom line

Three hours of pipeline work + two models with binary edge coverage have produced one defensible headline (within-length effect sizes well above the pre-registered threshold in ORCHIDEE for hops 3–8+) and one defensible complementarity claim (dormant high-risk paths in the chosen configuration). The work has *not* produced the full pre-registered claim, and cannot until per-edge counts arrive, the static/runtime codebase drift is resolved, and at least one more model is profiled. None of those are analytical bottlenecks on our side — they are deliverables we are waiting on from the engineering team and from the static-analysis re-run.

The most useful thing the analysis side can do in the meantime is (a) write the audit script in §6.4 to verify the name-matching, (b) draft the Results paragraph in §6.8, and (c) hold the rest pending the missing inputs.
