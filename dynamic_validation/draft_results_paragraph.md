# Draft text for the dynamic-profiling validation section

This draft is calibrated to the current state of the validation pipeline
(`summary.md`, this folder) and to the narrative position the placeholder
occupies in the manuscript: it sits in the Results section, immediately after
the per-model architectural-diversity paragraph (line 170) and immediately
before the "top 5% of functions ≈ 70–80% of risk" concentration result (line
172). Two paragraphs of Results, ~325 words. A matching Methods subsection
(~245 words) is drafted below. Word counts and tone follow the existing
paragraph in lines 152-161 (code churn) and the surrounding text.

> **Status (2026-05-17 refresh):** scoped to ORCHIDEE + PCR-GLOBWB.
> Per-edge runtime counts arrived in the latest export and Analyses 2, 3
> and 4 now produce numbers — folded into the two paragraphs below. VIC
> and HYPE runs are still pending; once they land, the marked passages
> below need refreshing. See the "Update checklist when VIC and HYPE
> arrive" section at the end of this file for what specifically changes.

---

## To replace the placeholder at line 171

**Heading (optional):** *Dynamic profiling preferentially exercises statically high-risk paths.*

<!-- [UPDATE WHEN VIC AND HYPE ARRIVE]: change "two of the four large models"
to "four large models", drop the "for which both reproducible execution and
runtime instrumentation were available" qualifier if all four are covered,
and replace "(ORCHIDEE and PCR-GLOBWB)" with the full list. -->
We complement these static results with a dynamic profiling validation on
two of the four large models for which both reproducible execution and
runtime instrumentation were available (ORCHIDEE and PCR-GLOBWB). For each
model we ran a single developer-shipped reference configuration (Methods)
and recorded which static call-graph edges were traversed at runtime, then
asked whether the statically identified high-risk paths were preferentially
exercised. Because the path-level risk score $P_k$ aggregates over the
nodes of a path under a saturating-OR formulation (Eq. 5), it is monotone
in path length, and an unstratified top-vs-bottom decile contrast conflates
static risk with hops. To dissolve this confound we computed deciles of
$P_k$ within each path-length bin and compared coverage.
<!-- [UPDATE WHEN VIC AND HYPE ARRIVE]: if within-length effect sizes for
VIC and HYPE also land at A > 0.55 in the long-path bins, change "In
ORCHIDEE" below to "Across the four models" (or "in three of the four
models", etc.) and report a per-model A range rather than a single
trajectory. If they diverge, name the discrepancy and refer to a per-model
table in SM. The pre-registered decision rule (A > 0.55 in >=3 of 4 models;
see pre_registration.md §6) only becomes testable here. -->
In ORCHIDEE,
statically high-risk paths exercised more edges at runtime than
length-matched low-risk paths across every length bin of three or more
hops: the common-language effect size $A$ rose monotonically from 0.78
(hops 3–4) through 0.93 (hops 5–7) to 1.00 (hops 8+), with Wilcoxon
$p < 10^{-6}$ in all three bins and well above the $A > 0.55$ threshold
adopted for the churn analysis (Table S2). For paths of one or two hops the
framework provided no discrimination — at length one, high-risk single-
function calls were less exercised than low-risk ones, consistent with the
framework's design as an integrator over multi-node execution chains rather
than a per-function score (Supplementary Materials, Fig. SX).

Among the paths that runtime did exercise, static $P_k$ orders runtime
intensity. The sum of per-edge call counts along each path
($E_k^{sum}$) correlates strongly with the median static risk score
across the uncertainty ensemble (Spearman $\rho = 0.79$, 95 % bootstrap
CI [0.77, 0.81] over 2,100 exercised paths in ORCHIDEE; Kendall
$\tau = 0.62$). The same conclusion holds when we hold execution
frequency fixed. Restricting to the runtime-heavy decile and comparing
high- versus low-$P_k$ paths within that subset, cyclomatic complexity
remains significantly higher in the high-$P_k$ stratum (Wilcoxon
$p = 0.004$, $A = 0.59$). Mahalanobis nearest-neighbour matching of
high- and low-$P_k$ paths on $E_k^{sum}$, hops and total statement
count (211 pairs in ORCHIDEE) produces an even larger gap in
cyclomatic complexity (median 77 versus 12 per function; $A = 0.81$),
indicating that static topology adds information beyond what runtime
profiling alone provides. PCR-GLOBWB shows the same qualitative
direction but the number of exercised paths is too small for either
the rank correlation (CI brackets zero) or the matched contrast (fewer
than five matching candidates per stratum) to support quantitative
claims (Supplementary Materials).

The same data substantiate the framework's intended role as a
forward-looking indicator of latent vulnerability. Within each length
bin of three or more hops in ORCHIDEE, no top-decile high-risk path
had every edge exercised in the reference configuration, and at
intermediate lengths (hops 3–4) 5 % of within-length top-decile paths
had no edge exercised at all.
<!-- [UPDATE WHEN VIC AND HYPE ARRIVE]: if VIC/HYPE generate enough
within-length top-decile paths to quantify dormancy, replace the
PCR-GLOBWB sentence with per-model numbers and note that three of the
four models produced quantitative dormancy estimates. -->
Dynamic profiling thus captures only a subset of the execution paths
the static framework identifies as systemically risky. Together, the
preferential exercise of high-$P_k$ paths, the persistent gap in
cyclomatic complexity after matching on execution frequency, and the
substantial dormant fraction among statically high-risk paths
corroborate the framework on three complementary axes: paths that are
exercised at runtime are exercised preferentially when their static
$P_k$ is high; static topology continues to discriminate after
controlling for execution intensity; and the high-risk paths that are
not exercised are precisely those whose surveillance the method is
designed to support.

---

## To add as a new Methods subsection (after "The risk metrics", before any subsequent subsection)

**Heading:** *Dynamic profiling validation.*

<!-- [UPDATE WHEN VIC AND HYPE ARRIVE]: change "two of the four large models"
to "four large models"; add VIC and HYPE to the list, and add their reference
configurations (Livneh CONUS single-basin for VIC; eWaterCycle hype-grpc4bmi
demo for HYPE - see pre_registration.md §7). -->
For two of the four large models with both reproducible execution and
available runtime instrumentation (ORCHIDEE and PCR-GLOBWB), we ran a
single developer-shipped reference configuration — a single-site offline
run with CRUJRA forcing for ORCHIDEE; the GMD 2018 5 arc-min global run
for PCR-GLOBWB (Supplementary Materials) — and recorded the set of
caller-to-callee edges traversed by the profiler. Configurations were
fixed independently of the static path rankings to guard against
circularity, following a pre-registered protocol locked before any
analysis was run against the dynamic data (`dynamic_validation/`,
analysis repository). A static call-graph edge was treated as
runtime-exercised if it appeared in the recorded trace; a static path was
treated as exercised under a looser criterion (any edge exercised) and a
stricter criterion (all edges exercised). We report a third per-path
metric — the fraction of edges of the path that were exercised — which is
insensitive to path length. Because $P_k$ is monotone in path length under
the saturating-OR aggregation (Eq. 5), we computed deciles of $P_k$
within each path-length bin (1, 2, 3–4, 5–7 and 8+ hops) and compared
coverage by within-bin two-sample Wilcoxon rank-sum tests, reporting the
common-language effect size $A$ (Vargha and Delaney) as in Table S2.
Per-edge runtime call counts were aggregated to a path-level intensity
$E_k^{sum}$ by summing along the static path edges; an originally
pre-registered $E_k^{min}$ (the bottleneck aggregator) collapsed to a
near-binary indicator on first contact with the data and was demoted to
a robustness check, reported alongside the Spearman correlation.
Rank correlation between $P_k$ and $E_k^{sum}$ was estimated by Spearman
$\rho$ with a 1,000-resample percentile bootstrap CI; the reverse
predictive check (Analysis 3) restricted to the top decile of paths by
$E_k^{sum}$ and compared high- versus low-$P_k$ strata on path-mean
node attributes; the matched-frequency contrast (Analysis 4) paired
high- and low-$P_k$ paths via Mahalanobis nearest-neighbour matching
with replacement on $E_k^{sum}$, hops and total statement count
(`MatchIt`). We treat the framework as predictive within a given length
bin when $A > 0.55$ and the Wilcoxon test rejects at $p < 0.05$, the
same threshold used in the churn analysis; the pre-registered global
decision rule additionally requires Spearman $\rho > 0$ with bootstrap
CI strictly above zero. The full analysis pipeline, the joined
static/dynamic path-level data and the locked pre-registration of
decisions are available in the `dynamic_validation/` directory of the
analysis repository.

---

## Supplementary materials and figures the draft refers to

These need to be added to the existing SM document before submission.

| Reference in draft | What to add | Source |
|---|---|---|
| "single-site offline ... GMD 2018 5 arc-min ..." | One Methods sub-section in SM listing the reference configurations per model, with sources (developer docs, ISIMIP protocol etc.) | `dynamic_validation/pre_registration.md` §7 already has this in table form. |
| "Fig. SX" (hops=1 reversal) | One supplementary figure: bar plot of frac_edge by hops bin × stratum × model | `dynamic_validation/results/synthesis/01b_within_length_deciles.pdf` |
| "with insufficient N ... (Supplementary Materials)" | Sentence noting the PCR-GLOBWB N limitation and the bins that were exercisable | `dynamic_validation/summary.md` §5.3 |
| "high-risk paths ... 5% had no edge exercised" | Supplementary table: within-length dormancy counts per model | `dynamic_validation/results/synthesis/05_dormant_summary.csv` |
| "Spearman $\rho = 0.79$ ..." | Supplementary figure: per-model scatter of $P_k$ vs $\log_{10}(E_k^{sum} + 1)$ | `dynamic_validation/results/synthesis/02_scatter_ORCHIDEE.pdf`, `02_scatter_PCR-GLOBWB.pdf`, table `02_rank_correlation.csv` |
| "Wilcoxon $p = 0.004$, $A = 0.59$" (reverse predictive) | Supplementary figure: boxplots of path-mean static metrics in high vs low $P_k$ within runtime-heavy decile | `dynamic_validation/results/synthesis/03_reverse_predictive_boxplots.pdf`, table `03_reverse_predictive.csv` |
| "211 pairs ... $A = 0.81$" (matched contrast) | Supplementary table: matched contrast statistics + SMD diagnostics | `dynamic_validation/results/synthesis/04_matched_contrast.csv`, `04_matched_smd_diagnostics.csv` |
| "dynamic_validation/ directory of the analysis repository" | Add a line in the data-availability section pointing to the directory | n/a |

---

## Things deliberately not claimed in the draft

So you can be sure they have not slipped in:

- That the pre-registered decision rule for "framework predicts at runtime"
  is met across all four target models. It is met *cleanly* for ORCHIDEE
  on both legs (Spearman with CI > 0; A > 0.55 in Analyses 1b, 3 and 4)
  but the "≥3 of 4 models" cardinality criterion cannot be evaluated until
  VIC and HYPE arrive. The draft scopes accordingly.
- That the codebase-drift question (Federico's email) is settled. The
  draft is silent on this; the assumption is that the engineering
  "phase 2" static re-run is being handled separately and either
  supersedes these numbers or confirms them. If phase 2 lands before
  submission, the numbers here may need refreshing.
- That PCR-GLOBWB contributes quantitatively to the validation. With only
  17 exercised paths, 5 runtime-heavy paths and 2 matching candidates per
  stratum it does not, and the draft says so explicitly.
- That VIC and HYPE have contributed. They have not, and the draft scopes
  honestly to "two of the four large models".

---

## Sentences to consider for the Discussion

Two short additions to the existing Discussion would tie the new results to
the framing already in lines 236-243:

After line 243 ("...audit the source code to assess transparency and call
chain vulnerabilities."), one sentence:

> *Empirically, the dynamic profiling validation in Fig. SX bears this out:
> high-risk paths that runtime profiling does exercise are exercised
> preferentially over length-matched low-risk paths, while a substantial
> fraction of statically high-risk paths are not fully traversed in any
> single reference configuration — precisely the latent vulnerabilities for
> which the static framework is designed to provide visibility.*

And, as a limitation in the limitations paragraph (around line 244-258),
one sentence:

> *The dynamic profiling validation reported here covers two configurations
> on two models; broader generalisation will require additional models and
> additional configurations per model, which we leave to follow-up work.*

<!-- [UPDATE WHEN VIC AND HYPE ARRIVE]: this limitation sentence either
disappears or shrinks to "additional configurations per model". -->

---

## Update checklist when VIC and HYPE arrive

When Federico's team ships runtime profiling for VIC and HYPE, the
following needs to happen, in order:

**1. Pipeline.** Drop the raw exports under `dynamic_validation/VIC/` and
`dynamic_validation/HYPE/`, then append two rows to `model_registry` in
`analysis/00_setup.R`:

```r
model_registry <- data.table(
  model      = c("ORCHIDEE", "PCR-GLOBWB", "VIC", "HYPE"),
  raw_subdir = c("ORCHIDEE", "PCR",        "VIC", "HYPE"),
  lang       = c("fortran",  "python",     "fortran", "fortran")
)
```

(Adjust `lang` for HYPE/VIC if their runtime traces are Python-side
instead.) Then re-run the pipeline:

```sh
Rscript dynamic_validation/analysis/01_join_static_dynamic.R
Rscript dynamic_validation/analysis/02_coverage.R
Rscript dynamic_validation/analysis/06_dormant_paths.R
Rscript dynamic_validation/analysis/07_cross_model_synthesis.R
```

The within-length deciles in `01b_within_length_deciles.csv` will now
include four models. Read the per-model A and Wilcoxon p from there.

**2. summary.md** (in this folder).

- Section 1 "Models with runtime data shipped so far" → update.
- Section 3 "Headline results" tables → re-populate from the refreshed
  `01b_within_length_deciles.csv` and `05_dormant_summary.csv`.
- Section 4 "How this builds on and defends the manuscript" → if all four
  models clear A > 0.55 in their long-path bins, replace "1 of 2 available
  models" wording with the actual count (e.g. "3 of 4 models"). The
  pre-registered decision rule becomes evaluable; state explicitly whether
  it is met.
- Section 5.3 (PCR-GLOBWB N limitation) and 5.4 (two models missing) →
  remove or trim.
- Section 6.3 (push for VIC and HYPE) → strike.

**3. draft_results_paragraph.md** (this file). Every comment block of the
form `<!-- [UPDATE WHEN VIC AND HYPE ARRIVE] -->` flags a sentence to
revise. There are five such blocks: two in the Results paragraphs, one in
the Methods subsection, one in the Discussion-limitations addition, and
this checklist itself. After updating, delete the comment markers so the
draft is ready to paste.

**4. pre-registered decision rule, the key call.** The rule
(`pre_registration.md` §6) requires A > 0.55 in ≥3 of 4 models AND
Spearman ρ > 0 (CI > 0). With four models in hand, the A criterion is
testable. Spearman remains blocked until the `runtime_calls` column lands
on the matched-static edge set; that block is independent of the
VIC/HYPE addition. If A is met but Spearman isn't, state the partial
support explicitly in the Results.

**5. Supplementary materials.** Refresh the figures listed in the
"Supplementary materials and figures the draft refers to" table above
from the regenerated outputs in `dynamic_validation/results/synthesis/`.
