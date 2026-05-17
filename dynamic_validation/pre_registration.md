# Pre-registration: dynamic validation of path-level software risk

**Paper:** Puy et al., *The topology of software risk in scientific models* (target: *Nature Computational Science*).
**Date locked:** 2026-05-16.
**Scope:** Dynamic profiling validation for four models — PCR-GLOBWB, VIC, ORCHIDEE, HYPE.

The decisions below are committed *before* any analysis is run against the dynamic data. They follow §6 of `handover_dynamic_validation.md`. Any deviation must be logged at the bottom of this file with date and reason.

---

## 1. Path-level aggregation of $E_k$

Let $e_{ij}$ be the runtime call count on directed edge $(i,j)$ from the profiler. For each static path $k$ comprising edges $\{(i,j)\}_k$:

- **Primary:** $E_k = \min_{(i,j) \in k}\, e_{ij}$ — the bottleneck-style "is the full path actually traversed".
- **Robustness:** $E_k^{\text{sum}} = \sum_{(i,j) \in k}\, e_{ij}$ — pipeline-style; more permissive.

Geometric mean is not used (handover lists it as a possible option; we exclude it to keep the report cleaner).

If any edge in the path is absent from the runtime edge set, its count is treated as $0$. Under the minimum aggregator this propagates to $E_k = 0$; under the sum aggregator it contributes nothing.

## 2. "Exercised" threshold for binary coverage (Analysis 1)

A path is *exercised* under two definitions, both reported:

- **Looser:** at least one of its edges appears in the runtime edge set.
- **Stricter:** every edge in the path appears in the runtime edge set.

(Note: the handover phrases this in terms of functions called rather than edges traversed. We operationalise as edges because the profiler output is edge-keyed; the two are equivalent when the static call graph faithfully represents the function-level structure.)

## 3. Stratum definition for high-risk / low-risk paths

- $P_k$ is the median of the per-path risk score across the $(\alpha, \beta, \gamma, p)$ ensemble produced by `softwareRisk::uncertainty_fun`.
- High-risk = top decile of $P_k$ **AND** $\text{CV}(P_k)$ below the model-specific median CV. This matches the high-risk filter used in the manuscript's bug-injection protocol for consistency.
- Low-risk = bottom decile of $P_k$ (no CV filter; the bottom decile is intrinsically lower-variance).
- Top/bottom deciles are computed within each model separately. No pooling across models.

Both risk forms in the analysis repo (additive $p=1$, compensatory $p=-1$) are reported. Primary results use the additive form; compensatory reported as robustness.

## 4. Effect size convention

Common-language effect size $A$ (Vargha & Delaney), matching Table S2 of the manuscript (churn analysis). For paired comparisons in Analysis 4 we report the paired-sample analogue.

## 5. Bootstrap procedure

Spearman $\rho$ and Kendall $\tau$ in Analysis 2: 1,000 path-level resamples (with replacement), percentile method 95% CI. Bootstrapping happens *within* the exercised subset per model — not across models.

## 6. Decision rule for "framework predicts at runtime"

The paper claims support if, jointly:

- Analysis 2: Spearman $\rho$ is positive with lower bound of bootstrap CI $> 0$.
- Analyses 3 & 4: common-language $A > 0.55$ for at least one of cyclomatic complexity, in-degree, or betweenness, in at least 3 of 4 models.

If one criterion is met and the other is not, the paper reports it as partial support and qualifies the claim accordingly. If neither is met, the paper reports the null and develops the dormant-paths framing (Analysis 5) as the substantive finding.

## 7. Input configurations (frozen — see handover §4)

| Model | Configuration | Source |
|---|---|---|
| PCR-GLOBWB | GMD 2018 paper global 5 arc-min run | Developer-published config dir |
| VIC | Livneh CONUS 1950–2013, single basin (Tuolumne or Upper Colorado) | Standard VIC paper input |
| ORCHIDEE | Single-site offline run (FLUXNET or SAFRAN), CRUJRA forcing | TRENDY protocol |
| HYPE | eWaterCycle containerised demo, `hype-grpc4bmi:feb2021` | Reproducible container, BMI interface |

Caveat to disclose in Methods: HYPE and ORCHIDEE configurations are smaller in scope than production setups (WWH; TRENDY global). This is a lower bound on runtime path coverage, not a ceiling.

## 8. What counts as a deviation

Any change to:
- The configurations in §7,
- The path-level aggregation in §1,
- The stratum definition in §3,
- The decision rule in §6,

must be recorded below with date and one-paragraph reason.

---

## Deviation log

**2026-05-17 — Primary E_k aggregator switched from `min` to `sum`.**
On the first runtime export with per-edge counts (ORCHIDEE, PCR-GLOBWB),
the `min` aggregator collapsed: only 38 of 3,152 ORCHIDEE paths and 0 of
101 PCR-GLOBWB paths had every edge exercised at least once. Under `min`
that means $E_k = 0$ for 99% of paths, which is informationally identical
to the binary `all_exercised` flag (Analysis 1 frac_all). The Spearman
correlation that Analysis 2 calls for, and the runtime-heavy decile cuts
that Analyses 3 and 4 use, are not interpretable on a near-degenerate
distribution. The `sum` aggregator (originally listed as robustness, §1)
preserves the intended "runtime-heavy" notion and is what Analyses 2–4
now use as primary. `min` is reported as robustness in Analysis 2 only.

This is a re-ranking of the two pre-registered aggregators, not the
introduction of a third. The original `min` motivation ("is the path
fully traversed at least once") is now covered by the binary
`all_exercised` flag in Analysis 1.
