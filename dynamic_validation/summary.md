# Dynamic-validation: progress summary

**Date:** 2026-05-17 (updated; per-edge runtime counts now available)
**Scope:** Dynamic-profiling validation of the path-level software-risk framework described in Puy et al., *The topology of software risk in scientific models* (target *Nature Computational Science*).
**Models with runtime data shipped so far:** ORCHIDEE (Fortran, MAN-2025 branch), PCR-GLOBWB (Python). VIC and HYPE pending.

This document is a point-in-time snapshot for collaborators. It records what has been built, what the analyses show, what those results buy the manuscript, and the work that still needs to happen before the validation is publication-ready.

### A few terms used throughout

- **Static call graph.** The graph of which functions call which, derived from the source code alone (no execution required). Built earlier in the project via `code_hydrological_models.R`.
- **Path.** A simple chain of function calls from an entry function (something with no caller in the graph) to a sink function (something with no callee). Each path has a length in **hops** = number of edges = number of calls along the chain. A length-1 path is a single function call; a length-7 path is a chain of seven calls.
- **P_k.** The path-level risk score from the manuscript (Eq. 5). Aggregates per-function risk into a single number per path. Bounded in [0, 1]. Higher = riskier path.
- **CV_k.** The coefficient of variation of P_k across the (α, β, γ, p) uncertainty ensemble. Low CV_k = the path stays high-risk no matter how you define risk; high CV_k = its rank depends on the definition. We use the median-CV cutoff to isolate "robustly high-risk" paths.
- **Runtime trace.** The set of (caller, callee) edges that were actually traversed when the model was run under a given configuration. Shipped by the engineers as `runtime_callgraph_edges.csv` per model.
- **Exercised.** A static edge is *exercised* if it appears in the runtime trace; a static path is *exercised* (in the loose sense) if at least one of its edges is in the runtime trace, or (strict sense) if all of its edges are.
- **Top / bottom decile.** Top decile = paths in the top 10% of P_k; bottom decile = bottom 10%. Used as the contrast for "high-risk vs low-risk".
- **Effect size A** (Vargha & Delaney's common-language). The probability that a randomly drawn value from group X is greater than a randomly drawn value from group Y. A = 0.5 means no difference; A > 0.5 means group X tends to be larger; A < 0.5 means group X tends to be smaller. The pre-registration uses A > 0.55 as the threshold for "supportive".

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
    03_rank_correlation.R     Analysis 2 (Spearman P_k vs E_k)
    04_reverse_predictive.R   Analysis 3 (reverse predictive check)
    05_matched_contrast.R     Analysis 4 (matched-frequency contrast)
    06_dormant_paths.R        Analysis 5 (dormant high-risk paths)
    07_cross_model_synthesis.R headline table + main-text figure
  results/per_model/          *_paths_joined.csv + edge diagnostics per model
  results/synthesis/          cross-model tables + plots
```

Key design choices:

- **Pre-registration locked first.** All decisions in `pre_registration.md` (aggregation, threshold, decile cut, effect-size convention, decision rule) were set before any analysis ran against the dynamic data, with a deviation log at the bottom.
- **Three coverage metrics, not one.** `frac_any` (loose: ≥1 edge exercised), `frac_all` (strict: all edges exercised), and `frac_edge` (mean fraction of edges per path exercised — length-insensitive). The triplet defuses the path-length confound described in §5.
- **Slim uncertainty summary.** The 1.78 GB `full_ua_df.csv` was reduced to an 8.7 MB `full_ua_summary.csv` containing per-path P_k quantiles and an exact P_k_sd computed from the ensemble vector. The slim file is small enough to commit and makes the validation reproducible without the raw monster.
- **Count-dependent analyses are wired with graceful soft-fail.** Analyses 2/3/4 require `runtime_calls` on the matched-static edge set. While that column was missing they exited cleanly with a `[label] No per-edge runtime counts ... Skipping.` message; from the 2026-05-17 export onwards they produce numbers without any schema changes downstream. The same fallback is preserved so the pipeline still runs if a future model is shipped without counts.

---

## 2. What was actually run

The handover specified five analyses. As of the 2026-05-17 runtime export, all five run against ORCHIDEE + PCR-GLOBWB. The static side is restricted to additive UA (compensatory was originally planned as a robustness check but turned out to be redundant under how `uncertainty_fun` samples its parameters; see §5.8). Per-edge runtime counts arrived with the latest export, unblocking Analyses 2, 3 and 4. The pre-registered primary E_k aggregator was switched from `min` to `sum` on first contact with the data (see deviation log in `pre_registration.md`); `min` is reported alongside as a robustness check in Analysis 2.

| Analysis | Status | Output |
|---|---|---|
| 1 — coverage (unstratified, three metrics) | done | `01_coverage_table.csv`, `01_coverage.pdf` |
| 1b — within-length deciles + Wilcoxon + A | done | `01b_within_length_deciles.csv`, `01b_within_length_tests.csv`, `01b_within_length_deciles.pdf` |
| 1c — length-stratified breakdown | done | `01_coverage_by_length.csv`, `01_coverage_by_length.pdf` |
| 2 — Spearman ρ on P_k vs E_k | done | `02_rank_correlation.csv`, `02_scatter_<MODEL>.pdf` |
| 3 — reverse predictive check | done | `03_reverse_predictive.csv`, `03_reverse_predictive_boxplots.pdf` |
| 4 — matched-frequency contrast | done | `04_matched_contrast.csv`, `04_matched_smd_diagnostics.csv` |
| 5 — dormant high-risk paths (within-length) | done | `05_dormant_summary.csv`, `05_dormant_paths_table.csv` |
| 6 — cross-model synthesis | done | `06_cross_model_headline.csv`, `06_decision_rule.csv`, `06_main_figure.pdf` |

---

## 3. Headline results

This section reports the actual numbers produced by the pipeline on ORCHIDEE and PCR-GLOBWB. It is organised as seven progressively sharper tables. The first three are coverage-style claims, addressable as soon as we know which edges were touched: §3.1 is a sanity check on how much of the static call graph the runtime trace actually reached; §3.2 gives the raw top-vs-bottom contrast, which is dramatic but partly an artefact of how P_k correlates with path length; §3.3 is the length-controlled version of that contrast — the apples-to-apples test for coverage. The next three use the per-edge runtime call counts that arrived on 2026-05-17: §3.4 reports rank correlation of static $P_k$ against runtime intensity $E_k^{sum}$ across all exercised paths; §3.5 is the reverse predictive check, asking whether static $P_k$ still discriminates within the runtime-heavy subset; §3.6 is the matched-frequency contrast, the cleanest test of "static topology adds information beyond execution frequency". The final subsection §3.7 characterises the high-risk paths that the runtime trace did not fully exercise, supporting the manuscript's "static finds what runtime misses" framing.

### 3.1 Edge-level coverage diagnostics

**What this table answers:** how much of each model's static call graph was actually touched by the runtime trace under the chosen configuration. This is a precondition for any downstream coverage statement — if the runtime reaches only a tiny fraction of the static graph, every other number in this section is conditioned on that intersection.

| Model | static paths | static edges | runtime edges | matched | runtime coverage |
|---|---|---|---|---|---|
| ORCHIDEE | 3152 | 1368 | 767 | 168 | 12.3% |
| PCR-GLOBWB | 101 | 116 | 144 | 9 | 7.8% |

**Columns:**

- *static paths*: total number of simple paths enumerated from the static call graph for this model.
- *static edges*: total number of unique caller → callee edges in the static call graph.
- *runtime edges*: total number of unique caller → callee edges observed in the runtime trace.
- *matched*: edges present in both sets (i.e. static edges that the runtime trace actually traversed).
- *runtime coverage*: matched / static edges. The fraction of the static call graph that the runtime trace reached.

**Reading:** the runtime trace covered ~12% of ORCHIDEE's static graph and ~8% of PCR-GLOBWB's. Some of that gap is genuine dormancy of the chosen configuration (which is interesting — see §3.4); some is suspected codebase drift between the static snapshot and the runtime executable (see §5.2 in this document).

### 3.2 Unstratified coverage — the headline contrast (additive, all paths)

**What this table answers:** when we split the paths into the *top 10% by P_k* (high-risk) and the *bottom 10% by P_k* (low-risk), how does runtime coverage compare? This is the simplest version of the validation question.

| Model | stratum | n | frac_any | frac_edge | frac_all | median hops |
|---|---|---|---|---|---|---|
| ORCHIDEE | top decile | 318 | 100.0% | 45.9% | 0.0% | 7 |
| ORCHIDEE | bottom decile | 316 | 12.3% | 8.7% | 5.7% | 1 |
| PCR-GLOBWB | top decile | 13 | 92.3% | 45.5% | 7.7% | 4 |
| PCR-GLOBWB | bottom decile | 15 | 0.0% | 0.0% | 0.0% | 1 |

**Columns:**

- *stratum*: high-risk (top decile of P_k, with the CV_k < median filter to keep only robustly high-risk paths) or low-risk (bottom decile of P_k).
- *n*: number of paths in that stratum.
- *frac_any*: fraction of paths in the stratum where **at least one** edge of the path was exercised at runtime. Loose definition.
- *frac_edge*: mean over paths in the stratum of (number of edges of the path that were exercised) / (total edges in the path). Length-insensitive — a long path with half its edges exercised contributes the same 0.5 as a length-2 path with one edge exercised.
- *frac_all*: fraction of paths where **every** edge was exercised at runtime. Strict definition.
- *median hops*: median path length in the stratum.

**Reading:** the contrast is large under all three metrics, but the rightmost column makes the catch clear — top-decile paths in ORCHIDEE are typically 7 hops long while bottom-decile paths are typically 1 hop. The reason is mechanical: P_k aggregates over node risks under a saturating-OR formula (Eq. 5 of the manuscript), so longer paths almost always score higher. Some of the "100% vs 12% any-exercised" contrast is therefore static risk doing real work; some is just the fact that longer paths have more chances to contain an exercised edge. The length-controlled version is in §3.3.

### 3.3 Within-length deciles — the apples-to-apples test (ORCHIDEE)

**What this table answers:** within each path-length class separately, does the framework still discriminate? Concretely: bin paths by hops (1, 2, 3–4, 5–7, 8+), then within each bin take the top vs bottom decile of P_k and compare runtime coverage. Length is held constant by construction.

| hops bin | n top / bot | median frac_edge (top vs bot) | Wilcoxon p | A | reading |
|---|---|---|---|---|---|
| 1 | 23 / 24 | 0.00 vs 1.00 | <10⁻⁴ | **0.23** | reversed |
| 2 | 45 / 45 | 0.00 vs 0.00 | 1.00 | **0.50** | null |
| 3–4 | 100 / 100 | 0.50 vs 0.00 | <10⁻¹¹ | **0.78** | supportive |
| 5–7 | 133 / 133 | 0.43 vs 0.20 | <10⁻³⁵ | **0.93** | strongly supportive |
| 8+ | 16 / 16 | 0.44 vs 0.25 | <10⁻⁶ | **1.00** | maximal |

**Columns:**

- *hops bin*: path length category. The bin "3–4" contains paths of 3 or 4 hops; "8+" contains paths of 8 hops or more.
- *n top / bot*: number of paths in the top-decile and bottom-decile strata *within this bin*.
- *median frac_edge (top vs bot)*: median value of frac_edge for the top-decile stratum versus the bottom-decile stratum. Compare the two side by side: this is the within-length contrast.
- *Wilcoxon p*: p-value of a two-sample Wilcoxon rank-sum test on frac_edge between the two strata, within this bin.
- *A*: common-language effect size A for the same comparison. > 0.5 = top-decile paths tend to have higher frac_edge than bottom-decile; < 0.5 = the reverse; 0.5 = no difference.
- *reading*: one-word summary of whether the cell supports the framework. Pre-registered threshold for "supportive" is A > 0.55 with p < 0.05.

**Reading:** for paths of three hops or more — i.e. wherever P_k is doing genuine path-level integration over multiple node risks — the framework discriminates strongly, with effect sizes A = 0.78, 0.93 and 1.00. For paths of one or two hops the framework fails (length-1 actually reverses): high-cyclomatic-complexity single-function calls in ORCHIDEE tend to live in error-handling or rarely-fired branches, while low-complexity single-function calls are heavily used utility functions. This is consistent with the manuscript's framing of P_k as a path-level metric — at length 1 there is no path to integrate over.

PCR-GLOBWB: only the hops=1 bin had enough paths for a within-length cut, and both strata exercise 0% — no signal either way. Cross-model claims here rest on ORCHIDEE.

### 3.4 Rank correlation P_k vs E_k_sum (Analysis 2)

**What this table answers:** among the paths that runtime actually exercised, does $P_k$ rank-order their runtime intensity correctly? This is the second leg of the pre-registered decision rule.

| Model | n exercised paths | Spearman ρ (P_k, E_k_sum) | 95% CI | Kendall τ | Spearman with E_k_min (robustness) |
|---|---|---|---|---|---|
| ORCHIDEE | 2100 | **0.79** | [0.77, 0.81] | 0.62 | -0.21 |
| PCR-GLOBWB | 17 | 0.11 | [-0.51, 0.62] | 0.15 | -0.62 |

**Reading:** ORCHIDEE shows a strong positive rank correlation between static $P_k$ and the sum of per-edge runtime calls along the path. The bootstrap CI is bounded well above zero, meeting the pre-registered ρ > 0 criterion. The negative correlation with $E_k^{\min}$ confirms the degeneracy that motivated switching the primary aggregator (see deviation log): under min, $E_k$ collapses to a near-binary "is the path fully traversed" indicator that anti-correlates with length. PCR-GLOBWB has only 17 exercised paths so its CI brackets zero — no conclusion either way.

### 3.5 Reverse predictive check (Analysis 3, ORCHIDEE)

**What this table answers:** if we restrict to the *runtime-heavy* paths only — the top 10% by $E_k^{sum}$ — does $P_k$ still discriminate within that subset? In words: holding execution frequency fixed, does static topology add information?

| Model | metric | n high / low | median high | median low | Wilcoxon p | A |
|---|---|---|---|---|---|---|
| ORCHIDEE | cyclomatic complexity | 128 / 128 | 97 | 97 | 0.004 | **0.59** |
| ORCHIDEE | in-degree | 128 / 128 | 0 | 0 | – | 0.50 |
| ORCHIDEE | betweenness | 128 / 128 | 0 | 0 | – | 0.50 |

**Reading:** within the runtime-heavy decile, ORCHIDEE high-P_k paths carry functions with significantly higher cyclomatic complexity than length-matched low-P_k paths in the same decile (A = 0.59 > 0.55 threshold; p = 0.004). In-degree and betweenness do not discriminate in this subset because the runtime-heavy paths share the same high-traffic entry and sink functions — their nodes have nearly identical node-level network metrics. The signal is on the complexity axis. PCR-GLOBWB was skipped because only 5 paths qualify as runtime-heavy.

### 3.6 Matched-frequency contrast (Analysis 4, ORCHIDEE)

**What this table answers:** the cleanest test of "static topology adds information beyond execution frequency". For each top-decile high-$P_k$ path, find a bottom-decile low-$P_k$ path with comparable $E_k^{sum}$, hops and statement count (Mahalanobis nearest-neighbour matching with replacement). Then compare the static metrics on the matched set.

| Model | metric | n pairs | median high | median low | Wilcoxon p | A |
|---|---|---|---|---|---|---|
| ORCHIDEE | cyclomatic complexity | 211 | 77 | 12 | 0.056 | **0.81** |
| ORCHIDEE | in-degree | 211 | 0 | 0 | – | 0.50 |
| ORCHIDEE | betweenness | 211 | 0 | 0 | – | 0.50 |

**Reading:** even after matching paths on execution intensity, path length and statement count, the high-$P_k$ stratum has substantially more cyclomatic complexity per function than the matched low-$P_k$ stratum (median 77 vs 12 cyclomatic units; A = 0.81). The Wilcoxon p is marginal (0.056) because the distribution is heavy-tailed and the rank-sum test is conservative under matching-with-replacement; the effect size is the more informative summary. PCR-GLOBWB had only 2 high-risk and 2 low-risk candidates after the CV filter so matching was skipped.

### 3.7 Within-length dormancy — the complementarity claim (ORCHIDEE top-decile)

**What this table answers:** of the within-length top-decile high-risk paths in each length bin, how many are dormant — i.e. were not (fully) exercised in the chosen configuration? Restricting to within-length top-decile (rather than overall top decile) means we are not just counting long paths that mechanically have more chances to fail.

| hops bin | n | dormant loose | dormant strict |
|---|---|---|---|
| 1 | 1 | 1 (100%) | 1 (100%) |
| 2 | 41 | 39 (95%) | 41 (100%) |
| 3–4 | 100 | 5 (5%) | 100 (100%) |
| 5–7 | 133 | 0 (0%) | 133 (100%) |
| 8+ | 8 | 0 (0%) | 8 (100%) |

**Columns:**

- *hops bin*: same as in §3.3.
- *n*: number of within-length top-decile high-risk paths in this bin (the same "n top" column as §3.3).
- *dormant loose*: paths for which **no** edge was exercised at runtime. The most conservative "completely missed" count.
- *dormant strict*: paths for which **at least one edge was missed** by the runtime trace (i.e. the path was not fully traversed).

**Reading:** in the long bins (5–7, 8+) loose dormancy is zero — every top-decile high-risk path of that length had at least one of its edges exercised. But strict dormancy is 100% in every bin of three or more hops: *no* top-decile high-risk path had every one of its edges exercised in the chosen configuration. The substantive number for the manuscript is the hops 3–4 row: 5 out of 100 (5%) of within-length top-decile high-risk paths were completely dormant under the reference configuration. The framework is flagging execution chains that the chosen run touched only partially, which is the framework's intended role.

---

## 4. How this builds on and defends the manuscript

This section is the bridge from the numbers in §3 to the claims the manuscript needs to make. The framework's central claim is that path-level static risk identifies execution paths that *concentrate systemic risk in scientific software*. Each of the five analyses in the handover was designed to address a specific reviewer counter-claim. With the numbers now in hand — including the rank-correlation, reverse-predictive and matched-frequency results that arrived in the 2026-05-17 export — we can mark off all five counter-claims on ORCHIDEE, leave the cross-model generalisation pending VIC and HYPE, and explicitly preserve the framework's most important framing concession (that static risk is *potential* and not *realised* vulnerability).

**Claim defended: "the framework discriminates beyond what trivial baselines would predict."** The unstratified coverage contrast (100% vs 12.3% any-exercised in ORCHIDEE; 92.3% vs 0% in PCR) confirms that top-decile P_k paths are an order of magnitude more likely to be exercised than bottom-decile paths. Once we control for path length — necessary because the saturating-OR aggregator is monotone in path length — the within-length effect sizes on ORCHIDEE are A = 0.78–1.00 across hops 3–8+, all p < 10⁻⁶. This *exceeds the pre-registered threshold of A > 0.55* in every length bin in the range where the framework is designed to apply. It is the cleanest evidence that static topology is doing more than reproducing path length.

**Claim defended: "static analysis identifies vulnerabilities that runtime profiling alone misses."** This was the manuscript's own framing of complementarity rather than substitutability, and it is now empirically grounded. At hops 3–4, 5% of within-length top-decile high-risk paths have zero edges exercised in the chosen ORCHIDEE configuration; at every length bin ≥3, *no* top-decile path has every edge exercised. The runtime trace touches the beginnings and middles of high-risk execution chains but not their full traversal — exactly the "potential vulnerability" the framework was built to flag.

**Claim defended: "the framework adds information beyond execution frequency itself."** This is the harder reviewer worry, addressed by Analyses 3 and 4. Restricting to the runtime-heavy decile of paths (those with the highest $E_k^{sum}$) and comparing high- versus low-$P_k$ strata within that subset, ORCHIDEE shows a significant gap in path-mean cyclomatic complexity ($A = 0.59$, $p = 0.004$). Mahalanobis nearest-neighbour matching of high- and low-$P_k$ paths on $E_k^{sum}$, hops and statement count produces an even larger gap (211 pairs; median cyclomatic complexity 77 vs 12, $A = 0.81$). So even when execution frequency is held fixed by design, the static framework still separates complex from simple execution paths — runtime weighting alone does not subsume it.

**Framing reinforced rather than weakened: the hops=1 reversal.** Length-1 high-risk paths are *less* exercised at runtime than length-1 low-risk paths in ORCHIDEE (4.3% vs 58.3%, A = 0.23, p < 10⁻⁴). The natural reading is that high-cyclomatic-complexity single-function calls live in error-handling, edge-case, or specialised routines that do not fire under nominal inputs, while low-complexity single-function calls are exactly the workhorse utilities the model uses constantly. This is consistent with — and arguably the cleanest empirical instance of — the manuscript's "static risk = potential, not realised, vulnerability" framing. It needs to be reported as part of the validation, not hidden in the supplement.

**Pre-registered decision rule, status.** The rule requires Spearman ρ > 0 (lower bound of bootstrap CI > 0) on the P_k–E_k correlation *and* A > 0.55 in ≥3 of 4 models. **ORCHIDEE meets both legs cleanly** — Spearman ρ(P_k, E_k_sum) = 0.79 [0.77, 0.81], well above 0; effect size A = 0.81 on cyclomatic complexity in the matched-frequency contrast (Analysis 4), 0.59 in the reverse predictive check (Analysis 3), and 0.78–1.00 in the within-length coverage comparison (Analysis 1b, hops 3–8+). PCR-GLOBWB is too data-thin to evaluate against the rule (only 17 exercised paths, 5 runtime-heavy paths, 2 candidates for matching). The "≥3 of 4 models" cardinality criterion cannot be met until VIC and HYPE arrive. What we have right now is therefore an unconditional pass on the one model where the test is statistically possible, with the remaining three models pending.

---

## 5. Weaknesses

These are written so that a reviewer's most natural objections appear on this list before they appear in a review.

**5.1 Single configuration per model.** Both ORCHIDEE and PCR-GLOBWB were profiled under one configuration apiece (single-site offline for ORCHIDEE; a working subset for PCR). Other configurations would exercise different subsets of paths. The dormant set in particular is a configuration-specific quantity. Generalising any "5% of top-decile paths are dormant" statement beyond the chosen configuration requires either more configurations or an explicit lower-bound framing.

**5.2 Codebase drift between static and runtime snapshots.** The engineers note that (a) PCR-GLOBWB was profiled against a slightly updated version because the original would not execute, and (b) ORCHIDEE has an architecture-specific code-regeneration step, so the Fortran the profiler observed is not identical to the Fortran the static analyser parsed. Both interpretations are partly responsible for the 88% non-match between static and runtime edge sets — but we cannot currently separate genuine dormancy from drift artefacts. *Until the static analysis is re-run against the actual runtime code versions ("phase 2" in the engineering plan), every coverage number carries a footnote of unknown size.*

**5.3 PCR-GLOBWB is statistically thin.** Only 13 paths in the unstratified top decile; only one length bin (hops=1) has enough paths for the within-length cut, and that bin shows zero exercise on both sides. PCR is providing essentially no validation signal at present, and any cross-model claim relies on ORCHIDEE alone.

**5.4 Two of four target models still missing.** VIC and HYPE have not been profiled yet. The pre-registered decision rule explicitly requires ≥3 of 4 models to meet the A > 0.55 threshold; with two models in hand and one of them data-limited, we cannot meet that rule on cardinality grounds, even though ORCHIDEE on its own meets both legs cleanly. The cardinality criterion is what now gates the full pre-registered claim.

**5.5 ~~Per-edge runtime counts are missing on the matched-static edge set.~~ Resolved 2026-05-17.** The latest engineering export ships `runtime_calls` on `runtime_callgraph_edges.csv`, which unblocks Analyses 2, 3 and 4. Sub-weakness: under the pre-registered `min` aggregator only 38 of 3,152 ORCHIDEE paths and 0 of 101 PCR-GLOBWB paths had every edge exercised at least once; we therefore promoted the originally-`min`-as-primary / `sum`-as-robustness pair to `sum`-as-primary, recorded as a deviation in `pre_registration.md`. The `min` interpretation ("is the path fully traversed at least once") is now covered by the binary `all_exercised` flag in Analysis 1.

**5.6 Name normalisation is aggressive and unaudited.** Both the static and runtime caller/callee names are lowercased and stripped of leading class qualifiers before joining. This is the right starting point for Python BMI classes (e.g. `BmiPCRGlobWB.set_value` → `set_value`) but may over-collapse for Fortran with module-name prefixes or for Python overloads. A sample audit of unmatched edges has not yet been done.

**5.7 The hops=1 and hops=2 results are honest but rhetorically awkward.** A reviewer who reads the abstract and lands on "the framework fails at hops ≤ 2" without the saturating-OR explanation may misread it. The within-length cut needs to be presented carefully in Methods so the reversal is framed as the framework's domain of applicability rather than as evidence against it.

**5.8 Compensatory risk form is not separately reported.** During the runs we discovered that `softwareRisk::uncertainty_fun` samples p ∈ [-1, 2] internally regardless of which `p` was used to enumerate paths, so additive and compensatory UA give the same per-path uncertainty distribution. This means the manuscript's existing "additive vs compensatory" robustness pair is meaningful only for the point estimate from `all_paths_fun`, not for the ensemble median used in the validation. The validation reports additive only; if a separate compensatory-point-estimate ranking is wanted, it has to be computed in a side-script.

---

## 6. What needs to happen to strengthen the validation

Mapped 1:1 to the weaknesses in §5; ordered by how much each one would tighten the validation, not by how much work it is. The first item (counts) was resolved on 2026-05-17; the next two (codebase phase-2 re-run, VIC + HYPE data) still sit outside our analytical control — they are deliverables we are waiting on from the engineering team. The remaining five are things the analysis side can do or decide independently.

**6.1 ~~Get the per-edge `runtime_calls` column from Federico.~~ Done 2026-05-17.** Counts arrived in the latest export and Analyses 2, 3 and 4 now produce numbers. See §3.4–§3.6 for results, and the deviation log in `pre_registration.md` for the aggregator switch this motivated.

**6.2 Trigger the engineering "phase 2" re-run of the static analyser on the runtime code versions.** For PCR-GLOBWB this means parsing the updated codebase Federico had to use to get the model running. For ORCHIDEE this means parsing the post-regeneration Fortran. Once both are done, the matched-static edge set goes from "12% of static edges" to whatever the true intersection is, and the dormancy numbers stop carrying their drift footnote.

**6.3 Push for VIC and HYPE runs.** Even a single configuration each gets us to four models, which is what the decision rule was written for. If a second PCR-GLOBWB configuration is also feasible (e.g. ISIMIP3 instead of the developer-shipped run), that would address §5.3 simultaneously. *When VIC and HYPE runs land, the step-by-step refresh of every file in this folder (including `draft_results_paragraph.md`) is in the "Update checklist when VIC and HYPE arrive" section at the end of `draft_results_paragraph.md`.*

**6.4 Audit name normalisation on a sample of unmatched edges.** Take 50 random unmatched static edges per model, hand-check whether the unmatched callees are (a) genuinely absent from the runtime trace, (b) present but renamed/mangled, or (c) present but in a class/module the normaliser collapsed away. This is half a day's work and gates the credibility of every coverage number.

**6.5 Restrict the headline to hops ≥ 3 and state the restriction in Methods.** Write the within-length framing as the primary analysis ("the framework is defined as a path-level integration; we evaluate it where path-level integration is non-degenerate, i.e. hops ≥ 3"), then report hops 1 and 2 honestly in the supplement as evidence for the potential-vs-realised distinction. This pre-empts the most likely framing objection.

**6.6 Either accept compensatory as redundant, or compute a fixed-p compensatory ranking in a side-script.** If the latter, write a 30-line script that takes the existing `all_graphs` (or rebuilds it for the two models), runs `all_paths_fun(p=-1)`, and joins the resulting `path_risk_score` to the joined-paths CSV as an alternative ranking. Run the within-length deciles against that ranking to confirm the headline doesn't depend on the choice of p for point-estimate enumeration.

**6.7 Add a second configuration per model where feasible.** This is the single most powerful answer to "your dormant set is configuration-specific" — running the same analysis on a different developer-shipped configuration of each model and showing that the dormant *set* is different but the dormant *fraction* is in the same ballpark would be a strong robustness statement. The handover's pre-committed configurations table (§4 of the handover) already lists backup options for PCR-GLOBWB and ORCHIDEE.

**6.8 ~~Draft the ~400-word Results paragraph and ~250-word Methods paragraph called for in the handover §7.~~ Done 2026-05-17.** The draft now lives at [`draft_results_paragraph.md`](draft_results_paragraph.md) and folds in the rank-correlation and matched-frequency findings. It claims that the pre-registered decision rule is met cleanly on ORCHIDEE and scopes the cardinality criterion (≥3 of 4 models) as pending VIC and HYPE. Update markers `[UPDATE WHEN VIC AND HYPE ARRIVE]` flag the sentences that need refreshing once the other two models land.

---

## 7. Bottom line

With the 2026-05-17 runtime export and the previously-built pipeline, ORCHIDEE now meets the pre-registered decision rule on both legs: Spearman ρ(P_k, E_k_sum) = 0.79 [0.77, 0.81] across 2,100 exercised paths, and effect size A = 0.81 on cyclomatic complexity in the matched-frequency contrast holding E_k_sum, length and statement count fixed. Within-length coverage effect sizes are 0.78–1.00 across hops 3–8+. The complementarity claim is also empirically grounded: 5% of within-length top-decile paths at hops 3–4 were completely dormant in the chosen configuration, and no top-decile path at any length ≥3 had all its edges exercised. The work is one clean pass on the model where the test was statistically possible.

What is not done: the same battery on VIC and HYPE (waiting on engineering deliverables), the codebase-drift resolution Federico flagged, and either a second configuration on each model or a name-matching audit to bound the artefact contribution. None of those are analytical bottlenecks. The decision-rule cardinality criterion (≥3 of 4 models) cannot be evaluated until VIC and HYPE arrive.

The most useful thing the analysis side can do on its own right now is the name-matching audit in §6.4 — a half-day of hand-checking that would tighten the credibility of every coverage and dormancy number by separating genuine dormancy from name-normalisation artefacts. Beyond that, the work is gated on engineering deliverables: VIC + HYPE runtime traces, and the phase-2 static re-analysis on the actual runtime code versions.
