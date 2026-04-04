[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17962642.svg)](https://doi.org/10.5281/zenodo.17962642)
# The topology of software risk in scientific models

[Arnald Puy](https://www.arnaldpuy.com/), Federico Díaz, Ulrike Proske, Olivia Richards,
Seth N. Linga, Samuel Flinders, Carmen Aguiló-Rivera, and Fernando G. Tinetti.

This study proposes a framework to identify risky paths in scientific models: sequences of
function calls whose potential failure can cascade into other parts of the software. Node
risk is scored as a weighted combination of cyclomatic complexity, in-degree and
betweenness centrality, with weights $\alpha$, $\beta$, and $\gamma$ respectively. Path
risk is then aggregated across nodes via a power-mean controlled by the exponent $p$. 
We illustrate the method with a synthetic example and apply it to fourteen global 
models spanning ecology, environmental science, hydrology and climate research.

## Abstract

*Software risk is studied in engineering and computer science but remains largely unexamined in scientific computing. Existing approaches typically assess risk at the level of individual functions or modules and are tied to a particular conception of risk, despite evidence of failures arising along execution paths and risk admitting multiple conceptualizations. Here we present a network-based framework that combines software-quality metrics with call graphs to identify execution paths concentrating systemic risk. The method treats the definition of risk as uncertain by comparing path rankings across a continuum space of risk definitions, thus accommodating the diverse priorities that shape risk assessment in scientific modelling. By illustrating its potential on fourteen models spanning ecology, environmental sciences, hydrology and climate change, we reveal that risk concentrates in a small subset of paths. Our approach provides a transparent basis for guiding testing, refactoring and verification in scientific software.*

## Dependencies

All analyses rely on the `softwareRisk` R package (available on [CRAN](https://cran.r-project.org/web/packages/softwareRisk/index.html)), which provides the core functions
`all_paths_fun()` and `uncertainty_fun()`. Additional R packages are loaded automatically
via `sensobol::load_packages()` at the top of each script. The key packages are:
`data.table`, `tidyverse`, `tidygraph`, `igraph`, `ggraph`, `sensobol`, `softwareRisk`,
`cowplot`, `openxlsx`, `scales`, and `benchmarkme`.

## Models

The fourteen models used in this study are:

| Model | Full name | Language | Domain |
|---|---|---|---|
| [CTSM](https://github.com/ESCOMP/CTSM) | Community Terrestrial Systems Model | Fortran + Python | Ecology, Environmental science, Hydrology, Climate research |
| [CWatM](https://github.com/iiasa/CWatM) | Community Water Model | Python | Hydrology, Environmental science |
| [DBH](https://hydro.iis.u-tokyo.ac.jp/DBH/index_files/Page394.htm) | Distributed Biosphere-Hydrological Model | Fortran | Hydrology, Ecology, Environmental science |
| [GR4J](https://github.com/EdgarEspitia/GR4J) | Génie Rural à 4 paramètres Journalier | Python | Hydrology |
| [H08](https://github.com/h08model/H08) | H08 Global Hydrological Model | Fortran | Hydrology, Environmental science |
| [HBV](https://github.com/johnrobertcraven/hbv_hydromodel) | Hydrologiska Byråns Vattenbalansavdelning | Python | Hydrology |
| HydroPy | Hydrological model | Python + Fortran | Hydrology, Climate research |
| [HYPE](https://sourceforge.net/projects/hype/files/) | Hydrological Predictions for the Environment | Fortran | Hydrology, Environmental science |
| [MHM](https://zenodo.org/records/8279545) | Mesoscale Hydrologic Model | Fortran + Python | Hydrology |
| [ORCHIDEE](https://forge.ipsl.jussieu.fr/orchidee/browser/branches/ORCHIDEE-MICT/tags/ORCHIDEE_MICT_8.4.1) | Organising Carbon and Hydrology In Dynamic Ecosystems | Fortran | Ecology, Environmental science, Hydrology, Climate research |
| [PCR-GLOBWB](https://github.com/UU-Hydro/PCR-GLOBWB_model) | PCRaster Global Water Balance model | Python | Hydrology, Environmental science |
| [SACRAMENTO](https://github.com/NOAA-OWP/sac-sma) | Sacramento Soil Moisture Accounting Model | Fortran | Hydrology |
| [SWAT](https://swatplus.gitbook.io/docs/source-code) | Soil and Water Assessment Tool | Fortran + Python | Hydrology, Ecology, Environmental science |
| [VIC](https://github.com/UW-Hydro/VIC) | Variable Infiltration Capacity model | Fortran + Python | Hydrology, Environmental science, Climate research |

## Repository structure

```
.
├── datasets/
│   ├── call_metrics/                  # Function calls per model and language
│   ├── cyclomatic_complexity_functions.csv  # Cyclomatic complexity per function
│   ├── descriptive_statistics/        # Lines of code, files, functions per model
│   ├── git_logs/                      # Git commit history per function
│   └── results_per_function/          # Function-level metrics
├── functions/                         # Custom R functions sourced by all scripts
│   ├── extract_sa_fun.R
│   ├── get_legend_fun.R
│   ├── to_tex_list_fun.R
│   ├── unnest_paths_tbl_fun.R
│   └── unnest_us_fun.R
├── code_synthetic_example.*           # Analysis 1: synthetic example
├── code_uncertainty_analysis.*        # Analysis 2: uncertainty & sensitivity
├── code_hydrological_models.*         # Analysis 3: domain models
├── code_scalability_test.*            # Analysis 4: scalability benchmark
└── README.md
```

Each analysis is available as `.R`, `.Rmd`, and `.pdf`.

## Analyses

### 1. Synthetic example (`code_synthetic_example`)

Constructs a small synthetic call graph to illustrate the method and verify its internal
consistency. Shows how path risk scores are computed and how the uncertainty and sensitivity
analysis works on a controlled example.

### 2. Uncertainty and sensitivity analysis (`code_uncertainty_analysis`)

Demonstrates the uncertainty and sensitivity analysis in depth, exploring how path risk
rankings respond to variation in the weights $\alpha$, $\beta$, $\gamma$ and the
aggregation exponent $p$.

### 3. Application to domain models (`code_hydrological_models`)

Applies the full pipeline to the fourteen models listed above. Includes:

- Construction of call graphs from raw call-metric datasets.
- Node-level risk scoring using cyclomatic complexity, in-degree, and betweenness
  centrality.
- Enumeration of all execution paths and their risk scores.
- Uncertainty and sensitivity analysis of path risk rankings.
- Analysis of metric robustness (swapping cyclomatic complexity for SLOC, betweenness for
  eigenvector centrality).
- Git-log analysis: whether high-risk functions exhibit higher churn.
- Real-model scalability timing and rank-stability assessment.

### 4. Scalability benchmark (`code_scalability_test`)

Stress-tests the method on synthetic layered DAGs of increasing size (up to ~2,500 nodes).
Reports wall-clock times for graph generation, path enumeration and uncertainty analysis,
as well as rank-stability metrics (top-1 stability, top-10 Jaccard similarity) across
graph sizes.

## Datasets

| Folder / file | Description |
|---|---|
| `call_metrics/` | CSV files with `from`/`to` function-call pairs, one file per model and language |
| `cyclomatic_complexity_functions.csv` | Cyclomatic complexity $C$ for every function in every model |
| `descriptive_statistics/` | Lines of code, number of files, functions, and modules per model |
| `results_per_function/` | File- and function-level metrics (LOC, bugs, complexity category) |
| `git_logs/` | Git commit history linked to individual functions, used for churn analysis |

## Functions

All custom functions live in the `functions/` folder and are sourced automatically at the
top of each script. Users do not need to load them manually.

| File | Purpose |
|---|---|
| `extract_sa_fun.R` | Extracts sensitivity indices from uncertainty analysis output |
| `get_legend_fun.R` | Utility to extract ggplot legends for composite figures |
| `to_tex_list_fun.R` | Formats path lists as LaTeX output |
| `unnest_paths_tbl_fun.R` | Unnests nested path tables from `all_graphs` |
| `unnest_us_fun.R` | Unnests uncertainty/sensitivity results from `all_graphs` |

## Replication

To reproduce all results:

1. Install `softwareRisk` and all packages listed in the Dependencies section.
2. Set your working directory to the root of this repository (the folder containing this
   README).
3. Run the scripts in order (1 → 2 → 3 → 4), or knit the corresponding `.Rmd`
   files. Each script sources all functions from `functions/` automatically.

> The code must be run from the repository root so that relative paths to `datasets/` and
> `functions/` resolve correctly.

## Citation

If you use this workflow, please cite:

A. Puy, F. Díaz, U. Proske, O. Richards, S. N. Linga, S. Flinders, C. Aguiló-Rivera, 
F. G. Tinetti (2026). Code and Datasets of The Topology of Software Risk in
Scientific Models. Zenodo. doi:10.5281/zenodo.17962642.

## License

MIT License

Copyright (c) 2026 Arnald Puy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
